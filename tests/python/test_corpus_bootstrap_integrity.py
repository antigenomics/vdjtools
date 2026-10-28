"""Published corpus inputs must remain complete and immutable across cache users."""
import hashlib
import json

import pytest

from vdjtools.signature import corpus as C


def _source(tmp_path, monkeypatch, payload=b'arrays'):
    res, source, cache = (tmp_path / name for name in ('res', 'source', 'cache'))
    res.mkdir(exist_ok=True)
    source.mkdir(exist_ok=True)
    monkeypatch.setenv('VDJTOOLS_CORPUS_DIR', str(cache))
    parts = {'npz': payload, 'json': b'{"format": 1}'}
    for ext, data in parts.items():
        (source / f'vsig_demo.{ext}').write_bytes(data)
    index = {'demo': {ext: {'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data)}
                      for ext, data in parts.items()}}
    (res / C.CORPORA_INDEX).write_text(json.dumps(index))
    return res, source


def test_resolver_repairs_partial_and_corrupt_cached_pair(tmp_path, monkeypatch):
    res, source = _source(tmp_path, monkeypatch)
    path = C.fetch_artifact('demo', sig='vsig', res_dir=res, repo='x/y',
                            base_url=source.as_uri(), quiet=True)
    download = C._download
    calls = []

    def local_download(url, target, digest, size):
        calls.append(target.suffix)
        return download((source / target.name).as_uri(), target, digest, size)

    monkeypatch.setattr(C, '_download', local_download)
    path.with_suffix('.json').unlink()
    assert C.resolve_artifact('demo', sig='vsig', res_dir=res, repo='x/y') == path
    assert path.with_suffix('.json').is_file()
    assert calls == ['.json']
    path.write_bytes(b'corrupt')
    C.resolve_artifact('demo', sig='vsig', res_dir=res, repo='x/y')
    assert path.read_bytes() == b'arrays'


def test_new_index_preserves_old_resource_paths(tmp_path, monkeypatch):
    res, source = _source(tmp_path, monkeypatch)
    old = C.fetch_artifact('demo', sig='vsig', res_dir=res, repo='x/y',
                           base_url=source.as_uri(), quiet=True)
    _source(tmp_path, monkeypatch, payload=b'new arrays')
    new = C.fetch_artifact('demo', sig='vsig', res_dir=res, repo='x/y',
                           base_url=source.as_uri(), quiet=True)
    assert new != old
    assert old.read_bytes() == b'arrays'
    assert new.read_bytes() == b'new arrays'


def test_concurrent_processes_download_each_part_once(tmp_path, monkeypatch):
    import subprocess
    import sys

    res, source = _source(tmp_path, monkeypatch)
    counter = tmp_path / 'downloads.txt'
    code = '''
import sys,time
from pathlib import Path
from vdjtools.signature import corpus as C
original = C._download

def counted(url, target, digest, size):
    with Path(sys.argv[3]).open('a') as log:
        log.write(target.suffix + '\\n')
    time.sleep(0.05)
    return original(url,target,digest,size)

C._download = counted
print(C.fetch_artifact('demo',sig='vsig',res_dir=Path(sys.argv[1]),repo='x/y',
                      base_url=Path(sys.argv[2]).as_uri(),quiet=True))
'''
    processes = [subprocess.Popen([sys.executable, '-c', code, str(res), str(source), str(counter)],
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                 for _ in range(2)]
    try:
        outputs = [p.communicate(timeout=30) for p in processes]
        assert all(p.returncode == 0 for p in processes), outputs
        assert outputs[0][0] == outputs[1][0]
        assert counter.read_text().splitlines() == ['.npz', '.json']
    finally:
        for p in processes:
            if p.poll() is None:
                p.kill()
                p.wait(timeout=5)


def test_complete_cache_needs_no_writer_lock(tmp_path, monkeypatch):
    import filelock

    res, source = _source(tmp_path, monkeypatch)
    path = C.fetch_artifact('demo', sig='vsig', res_dir=res, repo='x/y',
                            base_url=source.as_uri(), quiet=True)
    def forbidden(*args, **kwargs):
        raise AssertionError('complete cache attempted a write lock')
    monkeypatch.setattr(filelock, 'FileLock', forbidden)
    assert C.fetch_artifact('demo', sig='vsig', res_dir=res, repo='x/y', quiet=True) == path


def test_active_writer_times_out_and_killed_writer_releases_lock(tmp_path, monkeypatch):
    import subprocess
    import sys
    import time
    from filelock import Timeout

    res, source = _source(tmp_path, monkeypatch)
    path = C.fetch_artifact('demo', sig='vsig', res_dir=res, repo='x/y',
                            base_url=source.as_uri(), quiet=True)
    path.with_suffix('.json').unlink()
    ready = tmp_path / 'ready'
    code = '''
import sys,time
from pathlib import Path
from filelock import FileLock
with FileLock(sys.argv[1]):
    Path(sys.argv[2]).touch()
    time.sleep(30)
'''
    owner = subprocess.Popen([sys.executable, '-c', code, str(path.parent / '.fetch.lock'), str(ready)])
    try:
        deadline = time.monotonic() + 10
        while not ready.exists() and time.monotonic() < deadline:
            time.sleep(0.02)
        assert ready.exists(), 'writer did not acquire the lock'
        with pytest.raises(Timeout):
            C.fetch_artifact('demo', sig='vsig', res_dir=res, repo='x/y',
                             base_url=source.as_uri(), quiet=True, lock_timeout=0.05)
        owner.kill()
        owner.wait(timeout=5)
        assert C.fetch_artifact('demo', sig='vsig', res_dir=res, repo='x/y',
                                base_url=source.as_uri(), quiet=True, lock_timeout=1) == path
        assert path.with_suffix('.json').is_file()
        assert not list(path.parent.glob('*.part'))
    finally:
        if owner.poll() is None:
            owner.kill()
            owner.wait(timeout=5)


def test_verified_legacy_cache_is_reused_without_network(tmp_path, monkeypatch):
    import shutil

    res, source = _source(tmp_path, monkeypatch)
    cache = C.corpus_cache_dir()
    cache.mkdir()
    for path in source.iterdir():
        shutil.copyfile(path, cache / path.name)
        path.unlink()
    path = C.fetch_artifact('demo', sig='vsig', res_dir=res, repo='x/y',
                            base_url=source.as_uri(), quiet=True)
    assert path != cache / 'vsig_demo.npz'
    assert path.read_bytes() == b'arrays'
    assert (cache / 'vsig_demo.npz').read_bytes() == b'arrays'
