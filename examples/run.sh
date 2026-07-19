#!/usr/bin/env bash
# Launch the aging analysis notebook in marimo's interactive editor.
#   examples/run.sh              # open the editor
#   examples/run.sh --headless   # (or any marimo edit flag) forwarded through
#
# To serve it read-only instead (run-only app, no editing), use:
#   marimo run examples/aging.py
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
exec marimo edit "$here/aging.py" "$@"
