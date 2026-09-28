#!/usr/bin/env bash
# Copy reviewed direct-NN result artefacts into this chapter.
# No simulation data or trajectory files are copied. The results directory is
# supplied explicitly so this SI-facing repository contains no local paths.
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: bash $0 <reviewed-results-directory>" >&2
  exit 2
fi
SOURCE=$1
CHAPTER="$(cd "$(dirname "$0")/.." && pwd)"

for name in \
  nearest_neighbor_ionpair_tau_1e_all_forcefields_logy.png
do
  cp "$SOURCE/$name" "$CHAPTER/figures/$name"
done
cp "$SOURCE/nearest_neighbor_ionpair_tau_1e_all_forcefields.csv" "$CHAPTER/data/nearest_neighbor_tau_1e_all_forcefields.csv"
# The production seed table records filesystem provenance. That is useful
# internally but unsuitable for SI; retain only portable scientific columns.
cut -d, -f1-6 "$SOURCE/nearest_neighbor_ionpair_tau_1e_all_forcefields_by_seed.csv" \
  > "$CHAPTER/data/nearest_neighbor_tau_1e_all_forcefields_by_seed.csv"
