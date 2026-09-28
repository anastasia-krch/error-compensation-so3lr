#!/usr/bin/env bash
# Submit with: TEMP=300 SYSTEM_FILE=/path/to/start.in FORCEFIELD_FILE=/path/to/emimntf2_apple.ff sbatch npt.slurm
#SBATCH --job-name=applep_npt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=5-00:00:00
#SBATCH --output=applep_npt_%j.out

set -euo pipefail

: "${TEMP:?Set TEMP in kelvin, e.g. TEMP=300}"
: "${SYSTEM_FILE:?Set SYSTEM_FILE to the initial AMS system file}"
: "${FORCEFIELD_FILE:?Set FORCEFIELD_FILE to the APPLE&P force-field file}"
: "${AMSRC:?Set AMSRC to the AMS shell initialisation script}"

WORK_ROOT="${WORK_ROOT:-$PWD/runs}"
WORKDIR="${WORKDIR:-${WORK_ROOT}/${TEMP}K/npt}"
NSTEPS="${NSTEPS:-4000000}"
DT_FS="${DT_FS:-0.5}"

[[ -r "$SYSTEM_FILE" ]] || { echo "Missing SYSTEM_FILE: $SYSTEM_FILE" >&2; exit 1; }
[[ -r "$FORCEFIELD_FILE" ]] || { echo "Missing FORCEFIELD_FILE: $FORCEFIELD_FILE" >&2; exit 1; }
[[ -r "$AMSRC" ]] || { echo "Missing AMSRC: $AMSRC" >&2; exit 1; }

source "$AMSRC"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
unset AMS_SWITCH_LOGFILE_AND_STDOUT SCM_LOGFILE

mkdir -p "$WORKDIR"
cd "$WORKDIR"

AMS_JOBNAME=npt "$AMSBIN/ams" --delete-old-results <<EOF
MolecularDynamics
  Barostat
    ConstantVolume False
    Equal None
    Pressure 100000
    Scale XYZ
    Tau 1000
    Type Berendsen
  End
  BinLog
    DipoleMoment False
    PressureTensor False
    Time False
  End
  CalcPressure False
  Checkpoint
    Frequency 400000
  End
  InitialVelocities
    Temperature ${TEMP}
    Type Random
  End
  NSteps ${NSTEPS}
  Thermostat
    Tau 200.0
    Temperature ${TEMP}
    Type Berendsen
  End
  TimeStep ${DT_FS}
  Trajectory
    SamplingFreq 50
    WriteBonds False
    WriteCharges False
    WriteEngineGradients False
    WriteMolecules False
    WriteVelocities True
  End
End

Task MolecularDynamics

LoadSystem
  File ${SYSTEM_FILE}
End

Engine ForceField
  Type Apple&P
  ForceFieldFile ${FORCEFIELD_FILE}
EndEngine
EOF
