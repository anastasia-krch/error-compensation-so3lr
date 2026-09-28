#!/usr/bin/env bash
# Submit with: TEMP=300 RUN=1 START_XYZ=/path/to/extracted_npt_frame.xyz sbatch run_nvt.slurm
#SBATCH --job-name=so3lr_nvt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=so3lr_nvt_%j.out

set -euo pipefail

: "${TEMP:?Set TEMP in kelvin, e.g. TEMP=300}"
: "${RUN:?Set RUN, e.g. RUN=1}"
: "${START_XYZ:?Set START_XYZ to the initial configuration}"
SO3LR_BIN="${SO3LR_BIN:-so3lr}"
WORK_ROOT="${WORK_ROOT:-$PWD/runs}"
WORKDIR="${WORKDIR:-${WORK_ROOT}/${TEMP}K/run_${RUN}/nvt}"
MD_CYCLES="${MD_CYCLES:-10000}"
MD_STEPS="${MD_STEPS:-100}"
DT_FS="${DT_FS:-0.5}"
SEED="${SEED:-52}"

[[ -r "$START_XYZ" ]] || { echo "Missing START_XYZ: $START_XYZ" >&2; exit 1; }
command -v "$SO3LR_BIN" >/dev/null || { echo "Cannot find SO3LR_BIN: $SO3LR_BIN" >&2; exit 1; }

mkdir -p "$WORKDIR"
cd "$WORKDIR"

"$SO3LR_BIN" nvt \
  --input "$START_XYZ" \
  --output nvt.hdf5 \
  --log-file nvt.log \
  --temperature "$TEMP" \
  --dt "$DT_FS" \
  --md-cycles "$MD_CYCLES" \
  --md-steps "$MD_STEPS" \
  --nhc-chain 3 \
  --nhc-steps 2 \
  --nhc-thermo 200.0 \
  --seed "$SEED" \
  --precision float32 \
  --lr-cutoff 12.0 \
  --dispersion-damping 2.0 \
  --buffer-sr 1.25 \
  --buffer-lr 1.25 \
  --total-charge 0 \
  --save-buffer 50 \
  --no-relax \
  --restart-save nvt_restart.npz
