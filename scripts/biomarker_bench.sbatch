#!/bin/bash -l
#SBATCH --job-name=vdjt-biomarker
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Biomarker association / co-occurrence benchmark on Aldan-3 (appendix/bench_biomarker.py).
# One dataset per job:
#   aldan3 slurm submit scripts/biomarker_bench.sbatch --env vdjtools -- hip
#   aldan3 slurm submit scripts/biomarker_bench.sbatch --env vdjtools -- covid19
#   aldan3 slurm submit scripts/biomarker_bench.sbatch --env vdjtools -- covid19_vacc
# The env is a conda env with vdjtools installed from this checkout (pip install -e ".[dev,test]",
# which compiles vdjtools._core for the pgen model). Data read directly from the cluster.
set -euo pipefail
source /opt/conda/etc/profile.d/conda.sh 2>/dev/null || true
conda activate "${CONDA_ENV:-vdjtools}"
cd "${VDJT_DIR:-$HOME/vdjtools_bench}"
mkdir -p logs
PY="${PYTHON_BIN:-python}"
DATASET="${1:?usage: sbatch biomarker_bench.sbatch <hip|covid19|covid19_vacc>}"
FMBA=/projects/fmba_covid
VDJDB="${VDJDB:-/projects/immunestatus/vdjdb/vdjdb-2025-07-30/vdjdb.slim.txt}"
ORACLE="${ORACLE:-}"                                   # optional covid_associated_clonotypes.csv
echo "host=$(hostname)  job=${SLURM_JOB_ID:-local}  dataset=$DATASET  cpus=${SLURM_CPUS_PER_TASK:-?}"

case "$DATASET" in
  hip)   # local copy of Emerson HIP (counts + cmv/hla metadata); validate vs VDJdb-CMV
    "$PY" appendix/bench_biomarker.py --dataset hip --data-dir "$FMBA/hip_full" \
        --min-incidence 8 --vdjdb "$VDJDB"
    ;;
  covid19)   # 574 paired TRA+TRB subjects; COVID(current+past) vs healthy + α-β co-occurrence
    "$PY" appendix/bench_biomarker.py --dataset covid19 \
        --data-dir "$FMBA/data" --metadata "$FMBA/metadata_fmba_full.txt" --sample-col id \
        --min-incidence 10 ${ORACLE:+--oracle "$ORACLE"}
    ;;
  covid19_vacc)   # pre/post vaccination; sample id = 12-digit `name`
    "$PY" appendix/bench_biomarker.py --dataset covid19_vacc \
        --data-dir "$FMBA/vaccine/corr_func" \
        --metadata "$FMBA/vaccine/processed_metadata.tsv" --sample-col name --min-incidence 5
    ;;
  *) echo "unknown dataset: $DATASET" >&2; exit 2 ;;
esac
