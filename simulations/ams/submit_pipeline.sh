#!/usr/bin/env bash
# Usage: TEMP=300 RUN=1 NPT_RESTART_RKF=/path/to/MDStep4000000.rkf AMSRC=/path/to/amsbashrc.sh ./submit_pipeline.sh

set -euo pipefail

: "${TEMP:?Set TEMP in kelvin}"
: "${RUN:?Set RUN, e.g. RUN=1}"
: "${NPT_RESTART_RKF:?Set NPT_RESTART_RKF to the NPT checkpoint for this run}"
: "${AMSRC:?Set AMSRC}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_ROOT="${WORK_ROOT:-${SCRIPT_DIR}/runs}"
NVT_STEPS="${NVT_STEPS:-1000000}"
NVE_STEPS="${NVE_STEPS:-10000000}"
BASE="${WORK_ROOT}/${TEMP}K/run_${RUN}"

npt_engine="${NPT_ENGINE_RKF:-$(dirname "$NPT_RESTART_RKF")/ams.rkf}"
nvt_job=$(sbatch --parsable --export=ALL,TEMP="$TEMP",WORKDIR="$BASE/nvt",RESTART_RKF="$NPT_RESTART_RKF",ENGINE_RKF="$npt_engine",AMSRC="$AMSRC",NSTEPS="$NVT_STEPS" "$SCRIPT_DIR/nvt.slurm")
nve_restart="$BASE/nvt/nvt.results/MDStep${NVT_STEPS}.rkf"
nve_engine="$BASE/nvt/nvt.results/ams.rkf"
nve_job=$(sbatch --parsable --dependency="afterok:${nvt_job}" --export=ALL,TEMP="$TEMP",WORKDIR="$BASE/nve",RESTART_RKF="$nve_restart",ENGINE_RKF="$nve_engine",AMSRC="$AMSRC",NSTEPS="$NVE_STEPS" "$SCRIPT_DIR/nve.slurm")

printf 'NVT: %s\nNVE: %s\n' "$nvt_job" "$nve_job"
