#!/usr/bin/env bash
# vdjtools bootstrap — reproducible install.
#
# Steps:
#   1. Create/update the `vdjtools` conda environment (python + mmseqs2 + C++ toolchain).
#   2. (optional) editable-install co-developed sibling parents from ../ if present.
#   3. pip install -e ".[dev,test]" (compiles the _core C++ extension).
#
# Flags:
#   --no-conda       Skip conda env creation (use the already-active environment).
#   --dev-parents    Editable-install ../seqtree ../arda ../vdjmatch if they exist locally
#                    (they are early-alpha and co-developed; otherwise PyPI releases are used).
#   --tests          Run the fast test suites after install.
#
# Usage: bash setup.sh [--no-conda] [--dev-parents] [--tests]
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="vdjtools"
USE_CONDA=1
DEV_PARENTS=0
DO_TESTS=0

for arg in "$@"; do
  case "$arg" in
    --no-conda)    USE_CONDA=0 ;;
    --dev-parents) DEV_PARENTS=1 ;;
    --tests)       DO_TESTS=1 ;;
    *) echo "Unknown flag: $arg" >&2; exit 2 ;;
  esac
done

log() { printf '\033[1;34m[vdjtools]\033[0m %s\n' "$*"; }

# --- 1. conda environment --------------------------------------------------
if [[ "$USE_CONDA" -eq 1 ]]; then
  if ! command -v conda >/dev/null 2>&1; then
    echo "conda not found on PATH; install miniconda/anaconda or pass --no-conda." >&2
    exit 1
  fi
  if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    log "conda env '$ENV_NAME' exists — updating from environment.yml"
    conda env update -n "$ENV_NAME" -f "$ROOT/environment.yml" --prune
  else
    log "creating conda env '$ENV_NAME' from environment.yml"
    conda env create -f "$ROOT/environment.yml"
  fi
  PY="conda run -n $ENV_NAME python"
else
  PY="python"
fi

# --- 2. co-developed sibling parents (optional) ----------------------------
if [[ "$DEV_PARENTS" -eq 1 ]]; then
  for parent in seqtree arda vdjmatch; do
    if [[ -f "$ROOT/../$parent/pyproject.toml" ]]; then
      log "editable-install ../$parent"
      $PY -m pip install -e "$ROOT/../$parent"
    fi
  done
fi

# --- 3. editable install (builds the _core extension) ----------------------
log "pip install -e .[dev,test] (builds _core)"
$PY -m pip install -e "$ROOT[dev,test]"

# --- 4. verification -------------------------------------------------------
$PY -c "import vdjtools, vdjtools._core as c; print('vdjtools', vdjtools.__version__, '| _core', c.version())"

if [[ "$DO_TESTS" -eq 1 ]]; then
  log "running fast tests"
  $PY -m pytest "$ROOT/tests/python" -q
  log "C++ tests: cmake -S . -B build -DVDJTOOLS_TESTS=ON && cmake --build build && ctest --test-dir build"
fi

log "done."
