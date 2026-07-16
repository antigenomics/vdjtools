#!/bin/bash -l
#SBATCH --job-name=vdjt-biomarker
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Biomarker association / co-occurrence benchmark on Aldan-3 (appendix/bench_biomarker.py).
# Submit one dataset per job:
#   aldan3 slurm submit scripts/biomarker_bench.sbatch --env vdjtools -- hip
#   aldan3 slurm submit scripts/biomarker_bench.sbatch --env vdjtools -- covid19
#   aldan3 slurm submit scripts/biomarker_bench.sbatch --env vdjtools -- covid19_vacc
# The env is a conda env with vdjtools installed from this checkout (pip install -e ".[dev,test]",
# which compiles vdjtools._core for the pgen model). Data lives on the cluster (no HF needed for
# covid): /projects/fmba_covid (covid19 tables + metadata_fmba_full.txt) and the vaccine subdir.
set -euo pipefail
mkdir -p logs
PY="${PYTHON_BIN:-python}"
DATASET="${1:?usage: sbatch biomarker_bench.sbatch <hip|covid19|covid19_vacc>}"
FMBA=/projects/fmba_covid
VDJDB="${VDJDB:-/projects/immunestatus/vdjdb/vdjdb.slim.txt}"
echo "host=$(hostname)  job=${SLURM_JOB_ID:-local}  dataset=$DATASET  cpus=${SLURM_CPUS_PER_TASK:-?}"

case "$DATASET" in
  hip)
    # public HF dataset; metadata.txt carries cmv + hla
    "$PY" appendix/bench_biomarker.py --dataset hip --min-incidence 8 --vdjdb "$VDJDB"
    ;;
  covid19)
    "$PY" appendix/bench_biomarker.py --dataset covid19 \
        --data-dir "$FMBA/COV_V_usage_adjustment_v3/FMBA_functional" \
        --metadata "$FMBA/metadata_fmba_full.txt" --sample-col id \
        --oracle "$FMBA/covid_associated_clonotypes.csv" --min-incidence 10
    ;;
  covid19_vacc)
    "$PY" appendix/bench_biomarker.py --dataset covid19_vacc \
        --data-dir "$FMBA/vaccine/corr_func" \
        --metadata "$FMBA/vaccine/processed_metadata.tsv" --sample-col id --min-incidence 5
    ;;
  *) echo "unknown dataset: $DATASET" >&2; exit 2 ;;
esac
