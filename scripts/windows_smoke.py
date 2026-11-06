"""Exercise an installed Windows wheel: topology, native model, I/O and spawn."""
import tempfile
from pathlib import Path

import vdjtools
from vdjtools import io, stats
from vdjtools.cores import available_cores, physical_cores
from vdjtools.model import load_bundled
from vdjtools.model.generate import generate
from vdjtools.model.native import pgen_aa_batch


def summarize(frame):
    return stats.diversity_stats(frame)


def main():
    assert 1 <= physical_cores() <= available_cores()
    model = load_bundled("TRA", "isalgo")
    sample = generate(model, 8, seed=7, productive_only=True, engine="native", threads=0)
    scores = pgen_aa_batch(model, sample["junction_aa"].to_list(), threads=0)
    assert len(scores) == 8 and all(value > 0 for value in scores)
    with tempfile.TemporaryDirectory(prefix="vdjtools smoke ") as directory:
        path = Path(directory) / "sample.tsv"
        sample.write_csv(path, separator="\t", quote_style="never")
        results = io.map_samples(summarize, [("a", path), ("b", path)], workers=min(2, available_cores()), backend="process")
        assert len(results) == 2 and all(result[1].height == 1 for result in results)
    print("Installed wheel native model, I/O and spawned sample workers passed")


if __name__ == "__main__":
    main()
