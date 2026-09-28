#!/usr/bin/env bash
# Usage: TEMP=300 RUN=1 START_XYZ=/path/to/extracted_npt_frame.xyz ./submit_pipeline.sh

set -euo pipefail

: "${TEMP:?Set TEMP in kelvin}"
: "${RUN:?Set RUN, e.g. RUN=1}"
: "${START_XYZ:?Set START_XYZ}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_ROOT="${WORK_ROOT:-${SCRIPT_DIR}/runs}"
BASE="${WORK_ROOT}/${TEMP}K/run_${RUN}"
RUN_SEED="${SEED:-$((51 + RUN))}"

nvt_job=$(sbatch --parsable --export=ALL,TEMP="$TEMP",RUN="$RUN",SEED="$RUN_SEED",WORKDIR="$BASE/nvt",START_XYZ="$START_XYZ",SO3LR_BIN="${SO3LR_BIN:-so3lr}" "$SCRIPT_DIR/run_nvt.slurm")
nve_restart="$BASE/nvt/nvt_restart.npz"
nve_job=$(sbatch --parsable --dependency="afterok:${nvt_job}" --export=ALL,TEMP="$TEMP",RUN="$RUN",SEED="$RUN_SEED",WORKDIR="$BASE/nve",START_XYZ="$START_XYZ",RESTART_LOAD="$nve_restart",SO3LR_BIN="${SO3LR_BIN:-so3lr}" "$SCRIPT_DIR/run_nve.slurm")

printf 'NVT: %s\nNVE: %s\n' "$nvt_job" "$nve_job"
