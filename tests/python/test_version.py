"""The shipped version and the changelog's top entry must name the same release.

A release commit that writes the changelog and forgets `pyproject.toml` publishes under the previous
version, and a PyPI upload cannot be undone. Measured 2026-09-27, twice in one afternoon: arda 2.30.0
shipped its fix under the old literal and never reached PyPI, and mirpy's 4.0.0 commit left
`__version__` at 3.20.2. Nothing else in this suite looks at either number.
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_the_pyproject_version_matches_the_changelog_top_entry():
    # Read by regex rather than with `tomllib`, which is stdlib only from 3.11 while this package
    # supports 3.10 -- a test-only import of it fails *collection*, so the whole suite errors out
    # rather than one test failing. (Cost, measured 2026-09-28: exactly that, on both 3.10 CI jobs.)
    pyproject = Path(ROOT, "pyproject.toml").read_text()
    got = re.search(r'^version = "(\d+\.\d+\.\d+)"', pyproject, re.M)
    head = re.search(r"^## (\d+\.\d+\.\d+)", Path(ROOT, "CHANGELOG.md").read_text(), re.M)
    assert got, "pyproject.toml has no top-level `version = \"x.y.z\"`"
    assert head, "CHANGELOG.md has no `## <version>` entry"
    assert got.group(1) == head.group(1)


def test_the_installed_version_is_read_from_that_file_and_not_a_literal():
    """`vdjtools.__version__` comes from distribution metadata, deliberately -- a hand-copied literal
    drifted from the release version once already. This pins the mechanism, not the value: an
    editable install of a source tree reports whatever it was installed at, which is why the test
    above compares the file rather than the import."""
    src = Path(ROOT, "python/vdjtools/__init__.py").read_text()
    assert "_dist_version(\"vdjtools\")" in src
    assert not re.search(r'^__version__ = "\d', src, re.M)
