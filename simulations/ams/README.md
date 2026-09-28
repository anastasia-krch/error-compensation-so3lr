# APPLE&P / AMS setup

These Slurm templates implement the APPLE&P workflow used here:

1. NPT: 4,000,000 steps at 0.5 fs (2 ns), Berendsen thermostat and barostat.
2. NVT: 1,000,000 steps at 0.5 fs (500 ps), Nosé--Hoover-chain thermostat.
3. NVE: 10,000,000 steps at 0.5 fs (5 ns) for production analysis.

Run `npt.slurm` once per temperature. Its 400,000-step checkpoints are separated by 200 ps; use the five checkpoints from the second nanosecond as the starts of the five production branches. Set `AMSRC`, `SYSTEM_FILE`, and `FORCEFIELD_FILE` for the NPT run.

For example, from this directory:

```bash
export AMSRC=/path/to/amsbashrc.sh
export NPT_RESTART_RKF=/path/to/npt.results/MDStep4000000.rkf
TEMP=300 RUN=5 ./submit_pipeline.sh
```

`submit_pipeline.sh` submits the 500 ps NVT and 5 ns NVE stages for one selected checkpoint. `WORK_ROOT` defaults to `./runs`.
