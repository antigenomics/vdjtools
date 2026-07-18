#!/usr/bin/env bash
# vdjtools bootstrap — reproducible install into a repo-local .venv with uv.
#
# Portable: runs under bash OR zsh (bash setup.sh / zsh setup.sh / ./setup.sh).
#
# Steps:
#   1. Create/activate a repo-local .venv (uv if present, else python -m venv).
#   2. (optional) editable-install co-developed sibling parents from ../ if present.
#   3. pip install -e ".[dev,test]" — compiles the _core C++ extension via scikit-build-core.
#
# Flags:
#   --dev-parents  Editable-install ../seqtree ../arda ../vdjmatch if they exist locally
#                  (they are early-alpha and co-developed; otherwise the PyPI releases are used).
#   --tests        Run the fast test suites after install.
#
# Requirements: a C++ toolchain (Xcode Command Line Tools on macOS, build-essential on Linux)
# for the native _core extension; scikit-build-core + pybind11 + CMake are pulled in as build
# deps automatically. MMseqs2 (arda's aligner) is needed ONLY for the annotation / slow arda
# round-trip tests — `brew install mmseqs2` or `conda install -c bioconda mmseqs2`; it is not
# needed for Pgen, generation, EM inference, or any of the analytics.
#
# Usage: bash setup.sh [--dev-parents] [--tests]
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"   # $0, not ${BASH_SOURCE}: works in bash AND zsh
DEV_PARENTS=0
DO_TESTS=0

for arg in "$@"; do
  case "$arg" in
    --dev-parents) DEV_PARENTS=1 ;;
    --tests)       DO_TESTS=1 ;;
    --no-conda)    ;;  # accepted for backward-compat; conda is no longer used (no-op)
    *) echo "Unknown flag: $arg" >&2; exit 2 ;;
  esac
done

log() { printf '\033[1;34m[vdjtools]\033[0m %s\n' "$*"; }

# --- 1. repo-local .venv (uv preferred) ------------------------------------
VENV="$ROOT/.venv"
if command -v uv >/dev/null 2>&1; then
  PIP="uv pip"
  [ -d "$VENV" ] || { log "creating .venv with uv"; uv venv "$VENV"; }
else
  log "uv not found — using python -m venv + pip (install uv for faster installs: https://docs.astral.sh/uv/)"
  PIP="python -m pip"
  [ -d "$VENV" ] || { log "creating .venv"; python -m venv "$VENV"; }
fi
# shellcheck disable=SC1091
. "$VENV/bin/activate"   # activate script is bash/zsh compatible

# --- 2. co-developed sibling parents (optional) ----------------------------
if [ "$DEV_PARENTS" -eq 1 ]; then
  for parent in seqtree arda vdjmatch; do
    if [ -f "$ROOT/../$parent/pyproject.toml" ]; then
      log "editable-install ../$parent"
      $PIP install -e "$ROOT/../$parent"
    fi
  done
fi

# --- 3. editable install (builds the _core extension) ----------------------
log "$PIP install -e .[dev,test] (builds _core)"
$PIP install -e "$ROOT[dev,test]"

# --- 4. verification -------------------------------------------------------
python -c "import vdjtools, vdjtools._core as c; print('vdjtools', vdjtools.__version__, '| _core', c.version())"

if [ "$DO_TESTS" -eq 1 ]; then
  log "running fast tests"
  python -m pytest "$ROOT/tests/python" -q
  log "C++ tests: cmake -S . -B build -DVDJTOOLS_TESTS=ON && cmake --build build && ctest --test-dir build"
fi

log "done."
