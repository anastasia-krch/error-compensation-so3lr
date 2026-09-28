#!/usr/bin/env bash
#SBATCH --job-name=nn_reforming
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=rome
#SBATCH --time=12:00:00

# Usage: sbatch run_nearest_neighbour.sbatch <reduced-ionpair.xyz> <label> <out-dir>
set -euo pipefail
if [[ $# -ne 3 ]]; then
  echo "Usage: sbatch $0 <reduced-ionpair.xyz> <label> <out-dir>" >&2
  exit 2
fi

"${PYTHON:-python3}" \
  "$(dirname "$0")/compute_nearest_neighbour.py" \
  --xyz "$1" --label "$2" --out-dir "$3"
