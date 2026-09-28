#!/usr/bin/env bash
# Manifest-driven AMS single point or geometry optimisation at a fixed level of theory.
#
# Produces the AMS output that compute_pes_scan.py and compute_neutral_control.py read.
# Parameterised by environment variable rather than tied to a cluster or storage location,
# like the templates in simulations/.
#
#   MANIFEST=scan.tsv ENGINE=adf    XC=PBE0 BASIS=TZ2P ./run_ams_singlepoint.sh 0
#   MANIFEST=scan.tsv ENGINE=appleP APPLE_FF=/path/emimntf2_apple.ff ./run_ams_singlepoint.sh 0
#   MANIFEST=scan.tsv ENGINE=so3lr  SO3LR_UV_PROJECT=/path/to/uv/project ./run_ams_singlepoint.sh 0
#   MANIFEST=ref.tsv  ENGINE=adf  TASK=GeometryOptimization ./run_ams_singlepoint.sh 0
#
# Under Slurm the row index comes from SLURM_ARRAY_TASK_ID; otherwise pass it as $1.
#
# MANIFEST is a TSV with a header line, then one row per geometry:
#     relative_or_absolute_xyz <TAB> total_charge <TAB> tag
#
# The `tag=... xc=... basis=...` line echoed below is not decoration: charge_common.py's
# log index keys on it, and skips any log without it rather than guessing which geometry
# and level of theory a result belongs to.
set -euo pipefail

: "${MANIFEST:?MANIFEST must be set (TSV: xyz, charge, tag)}"
: "${OUTDIR:=$PWD/ams_runs}"
: "${LOGDIR:=$PWD/ams_logs}"

ENGINE="${ENGINE:-adf}"
TASK="${TASK:-SinglePoint}"
XC="${XC:-PBE}"
BASIS="${BASIS:-TZP}"
AMSBIN="${AMSBIN:-${AMSHOME:-}/bin}"
AMS_NPROC="${AMS_NPROC:-1}"

command -v "$AMSBIN/ams" >/dev/null 2>&1 || {
    echo "Cannot find the AMS driver at $AMSBIN/ams." >&2
    echo "Source your AMS environment first, or set AMSBIN." >&2
    exit 1
}

ROW_INDEX="${SLURM_ARRAY_TASK_ID:-${1:-0}}"
LINE=$(( ROW_INDEX + 2 ))                  # +2 skips the header
read -r XYZ CHARGE TAG < <(sed -n "${LINE}p" "$MANIFEST")
[[ -n "${TAG:-}" ]] || { echo "no manifest row at line $LINE of $MANIFEST" >&2; exit 1; }
if [[ "$XYZ" != /* ]]; then XYZ="$(dirname "$(readlink -f "$MANIFEST")")/$XYZ"; fi
[[ -f "$XYZ" ]] || { echo "missing geometry $XYZ" >&2; exit 1; }

if [[ "$ENGINE" == "appleP" ]]; then
    : "${APPLE_FF:?APPLE_FF must point at the APPLE&P force-field file}"
    [[ -f "$APPLE_FF" ]] || { echo "missing force field $APPLE_FF" >&2; exit 1; }
    LEVEL_XC="APPLE&P"; LEVEL_BASIS="ff"
    # APPLE&P has no automatic atom typing in AMS, so the geometry file must already
    # carry per-atom ForceField.Type entries.
    ENGINE_BLOCK="Engine ForceField
    Type APPLE&P
    ForceFieldFile $APPLE_FF
EndEngine"
elif [[ "$ENGINE" == "so3lr" ]]; then
    : "${SO3LR_UV_PROJECT:?SO3LR_UV_PROJECT must point at the SO3LR AMS uv project}"
    [[ -d "$SO3LR_UV_PROJECT" ]] || {
        echo "missing SO3LR uv project $SO3LR_UV_PROJECT" >&2; exit 1; }
    SO3LR_CALCULATOR="${SO3LR_CALCULATOR:-$(dirname "$(readlink -f "$0")")/so3lr_ams_calculator.py}"
    [[ -f "$SO3LR_CALCULATOR" ]] || {
        echo "missing SO3LR AMS calculator $SO3LR_CALCULATOR" >&2; exit 1; }
    LEVEL_XC="SO3LR"; LEVEL_BASIS="ASE"
    # The adapter exposes SO3LR's native partial_charges via AMS's ASE pipe, rather
    # than substituting an AMS population analysis.
    ENGINE_BLOCK="Engine ASE
    Type File
    File $SO3LR_CALCULATOR
    Arguments
        lr_cutoff = ${SO3LR_LR_CUTOFF:-12.0}
        dtype = '${SO3LR_DTYPE:-float64}'
    End
    Python
        Type Uv
        Uv $SO3LR_UV_PROJECT
    End
EndEngine"
else
    LEVEL_XC="$XC"; LEVEL_BASIS="$BASIS"
    if [[ "$XC" == PBE ]]; then XC_BLOCK="        GGA $XC"; else XC_BLOCK="        Hybrid $XC"; fi
    # CM5 is requested alongside the defaults: it is the only one of ADF's charge
    # analyses that behaves for both the imidazolium hydrogens and the CF3 groups.
    ENGINE_BLOCK="Engine ADF
    Basis
        Type $BASIS
        Core None
    End
    XC
$XC_BLOCK
        Dispersion Grimme3 BJDAMP
    End
    CM5 Yes
    NumericalQuality Good
    Symmetry NOSYM
EndEngine"
fi

WORKDIR="$OUTDIR/${LEVEL_XC}_${LEVEL_BASIS}/${TAG}"
rm -rf "$WORKDIR"; mkdir -p "$WORKDIR" "$LOGDIR"
LOG="$LOGDIR/ams_${TAG}_${LEVEL_XC}_${LEVEL_BASIS}.out"
XYZ_ABS=$(readlink -f "$XYZ")
cd "$WORKDIR"

{
    echo "tag=$TAG xyz=$XYZ_ABS charge=$CHARGE engine=$ENGINE xc=$LEVEL_XC basis=$LEVEL_BASIS task=$TASK"
    "$AMSBIN/ams" -n "$AMS_NPROC" <<EOR
System
    GeometryFile $XYZ_ABS
    Charge $CHARGE
End

Task $TASK

Properties
    Dipolemoment Yes
    Charges Yes
End

$ENGINE_BLOCK
EOR
} > "$LOG" 2>&1

# AMS writes its progress log inside ams.results/, not the working directory.
if grep -q "NORMAL TERMINATION" ams.results/ams.log; then
    # Engine ASE owns high-precision values in ase.rkf; stdout rounds charges to
    # four decimals.  Keep the RKF next to its provenance-bearing log so collectors
    # can always use the exact AMS result without reconstructing a work-directory path.
    if [[ -f ams.results/ase.rkf ]]; then
        cp ams.results/ase.rkf "$LOGDIR/ams_${TAG}_${LEVEL_XC}_${LEVEL_BASIS}.rkf"
    fi
    echo "OK: $TAG ($LEVEL_XC/$LEVEL_BASIS) -> $LOG"
else
    echo "FAILED: $TAG ($LEVEL_XC/$LEVEL_BASIS); see $LOG" >&2
    exit 1
fi
