#!/usr/bin/env bash
# Submit with: TEMP=300 RESTART_RKF=/path/to/npt.results/MDStep4000000.rkf sbatch nvt.slurm
#SBATCH --job-name=applep_nvt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=5-00:00:00
#SBATCH --output=applep_nvt_%j.out

set -euo pipefail

: "${TEMP:?Set TEMP in kelvin, e.g. TEMP=300}"
: "${RESTART_RKF:?Set RESTART_RKF to the selected NPT checkpoint}"
: "${AMSRC:?Set AMSRC to the AMS shell initialisation script}"

WORK_ROOT="${WORK_ROOT:-$PWD/runs}"
WORKDIR="${WORKDIR:-${WORK_ROOT}/${TEMP}K/nvt}"
ENGINE_RKF="${ENGINE_RKF:-$(dirname "$RESTART_RKF")/ams.rkf}"
NSTEPS="${NSTEPS:-1000000}"
DT_FS="${DT_FS:-0.5}"

[[ -r "$RESTART_RKF" ]] || { echo "Missing RESTART_RKF: $RESTART_RKF" >&2; exit 1; }
[[ -r "$ENGINE_RKF" ]] || { echo "Missing ENGINE_RKF: $ENGINE_RKF" >&2; exit 1; }
[[ -r "$AMSRC" ]] || { echo "Missing AMSRC: $AMSRC" >&2; exit 1; }

source "$AMSRC"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
unset AMS_SWITCH_LOGFILE_AND_STDOUT SCM_LOGFILE

mkdir -p "$WORKDIR"
cd "$WORKDIR"

AMS_JOBNAME=nvt "$AMSBIN/ams" --delete-old-results <<EOF
MolecularDynamics
  BinLog
    DipoleMoment False
    PressureTensor False
    Time False
  End
  CalcPressure False
  Checkpoint
    Frequency ${NSTEPS}
  End
  InitialVelocities
    Type FromFile
    File ${RESTART_RKF}
  End
  NSteps ${NSTEPS}
  Thermostat
    Tau 200.0
    Temperature ${TEMP}
    Type NHC
  End
  TimeStep ${DT_FS}
  Trajectory
    SamplingFreq 10
    WriteBonds False
    WriteCharges False
    WriteEngineGradients False
    WriteMolecules False
    WriteVelocities True
  End
End

Task MolecularDynamics

LoadSystem
  File ${RESTART_RKF}
End

LoadEngine ${ENGINE_RKF}
EOF
