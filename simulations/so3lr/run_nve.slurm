#!/usr/bin/env bash
# Submit with: TEMP=300 RUN=1 START_XYZ=/path/to/start.xyz RESTART_LOAD=/path/to/nvt_restart.npz sbatch run_nve.slurm
#SBATCH --job-name=so3lr_nve
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=5-00:00:00
#SBATCH --output=so3lr_nve_%j.out

set -euo pipefail

: "${TEMP:?Set TEMP in kelvin, e.g. TEMP=300}"
: "${RUN:?Set RUN, e.g. RUN=1}"
: "${START_XYZ:?Set START_XYZ to the initial configuration}"
: "${RESTART_LOAD:?Set RESTART_LOAD to the NVT restart}"
SO3LR_BIN="${SO3LR_BIN:-so3lr}"
WORK_ROOT="${WORK_ROOT:-$PWD/runs}"
WORKDIR="${WORKDIR:-${WORK_ROOT}/${TEMP}K/run_${RUN}/nve}"
MD_CYCLES="${MD_CYCLES:-100000}"
MD_STEPS="${MD_STEPS:-100}"
DT_FS="${DT_FS:-0.5}"
SEED="${SEED:-52}"

[[ -r "$START_XYZ" ]] || { echo "Missing START_XYZ: $START_XYZ" >&2; exit 1; }
[[ -r "$RESTART_LOAD" ]] || { echo "Missing RESTART_LOAD: $RESTART_LOAD" >&2; exit 1; }
command -v "$SO3LR_BIN" >/dev/null || { echo "Cannot find SO3LR_BIN: $SO3LR_BIN" >&2; exit 1; }

mkdir -p "$WORKDIR"
cd "$WORKDIR"

"$SO3LR_BIN" nve \
  --input "$START_XYZ" \
  --output nve.hdf5 \
  --log-file nve.log \
  --temperature "$TEMP" \
  --dt "$DT_FS" \
  --md-cycles "$MD_CYCLES" \
  --md-steps "$MD_STEPS" \
  --seed "$SEED" \
  --precision float32 \
  --lr-cutoff 12.0 \
  --dispersion-damping 2.0 \
  --buffer-sr 1.25 \
  --buffer-lr 1.25 \
  --total-charge 0 \
  --save-buffer 50 \
  --no-relax \
  --restart-load "$RESTART_LOAD" \
  --restart-save nve_restart.npz
