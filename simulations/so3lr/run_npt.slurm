#!/usr/bin/env bash
# Submit with: TEMP=300 START_XYZ=/path/to/start.xyz sbatch run_npt.slurm
#SBATCH --job-name=so3lr_npt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=4-00:00:00
#SBATCH --output=so3lr_npt_%j.out

set -euo pipefail

: "${TEMP:?Set TEMP in kelvin, e.g. TEMP=300}"
: "${START_XYZ:?Set START_XYZ to the initial configuration}"
SO3LR_BIN="${SO3LR_BIN:-so3lr}"
WORK_ROOT="${WORK_ROOT:-$PWD/runs}"
WORKDIR="${WORKDIR:-${WORK_ROOT}/${TEMP}K/npt}"
MD_CYCLES="${MD_CYCLES:-40000}"
MD_STEPS="${MD_STEPS:-100}"
DT_FS="${DT_FS:-0.5}"
SEED="${SEED:-52}"

[[ -r "$START_XYZ" ]] || { echo "Missing START_XYZ: $START_XYZ" >&2; exit 1; }
command -v "$SO3LR_BIN" >/dev/null || { echo "Cannot find SO3LR_BIN: $SO3LR_BIN" >&2; exit 1; }

mkdir -p "$WORKDIR"
cd "$WORKDIR"

"$SO3LR_BIN" npt \
  --input "$START_XYZ" \
  --output npt.hdf5 \
  --log-file npt.log \
  --temperature "$TEMP" \
  --pressure 1.0 \
  --dt "$DT_FS" \
  --md-cycles "$MD_CYCLES" \
  --md-steps "$MD_STEPS" \
  --nhc-chain 3 \
  --nhc-steps 2 \
  --nhc-thermo 200.0 \
  --nhc-baro 1000.0 \
  --seed "$SEED" \
  --precision float32 \
  --lr-cutoff 12.0 \
  --dispersion-damping 2.0 \
  --buffer-sr 1.25 \
  --buffer-lr 1.25 \
  --total-charge 0 \
  --save-buffer 50 \
  --no-relax \
  --restart-save npt_restart.npz
