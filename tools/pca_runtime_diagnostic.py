"""Separate exact PCA input/runtime failures from preceding cohort feature work."""
import faulthandler
import json
from pathlib import Path
import subprocess
import sys
import tempfile

faulthandler.enable()
import vdjtools  # initialize the supported native ceilings before numerical imports
import numpy as np
from vdjtools.signature import corpus
from threadpoolctl import threadpool_info


def report(stage, matrix):
    print(json.dumps({'stage': stage, 'shape': matrix.shape,
                      'finite': bool(np.isfinite(matrix).all()),
                      'pools': threadpool_info()}), flush=True)


if len(sys.argv) > 1:
    matrix = np.load(sys.argv[1])
    from sklearn import decomposition  # load the same pools as TablePCA
    report('fresh process exact input', matrix)
    corpus.fit_pca(matrix, 2, seed=9, copy=False)
    print('FRESH_PCA_PASS', flush=True)
else:
    import runpy
    original = corpus.fit_pca

    def inspect_fit(matrix, *args, **kwargs):
        report('parent exact input', matrix)
        with tempfile.TemporaryDirectory(prefix='vdjtools-pca-') as directory:
            path = Path(directory) / 'matrix.npy'
            np.save(path, matrix)
            result = subprocess.run([sys.executable, __file__, str(path)], timeout=120)
            print(f'fresh PCA exit: {result.returncode}', flush=True)
        report('parent before PCA', matrix)
        return original(matrix, *args, **kwargs)

    corpus.fit_pca = inspect_fit
    tests = runpy.run_path('tests/python/test_bounded_signature_fit.py')
    tests['test_sparse_table_pca_matches_dense_standard_scaler']()
    print('COHORT_PCA_PASS', flush=True)
