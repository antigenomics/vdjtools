"""Set signature kernel limits before importing the command implementation."""
import os
import sys


def main():
    if sys.argv[1:2] == ['signature']:
        from .cores import _KERNEL_ENV

        os.environ.update(dict.fromkeys(_KERNEL_ENV, '1'))
    from .cli import app

    app()
