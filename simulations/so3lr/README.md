# SO3LR setup

These Slurm templates reproduce the SO3LR sequence used for the simulations:

1. NPT: 4,000,000 steps at 0.5 fs.
2. NVT: 1,000,000 steps at 0.5 fs (500 ps) for each independent run.
3. NVE: 10,000,000 steps at 0.5 fs for production analysis.

Run `run_npt.slurm` once per temperature, then extract five configurations at 200 ps intervals from its second nanosecond. Set `START_XYZ` to one of these configurations and run the NVT--NVE pipeline:

```bash
export SO3LR_BIN=/path/to/so3lr
export START_XYZ=/path/to/start.xyz
TEMP=300 RUN=1 ./submit_pipeline.sh
```

`submit_pipeline.sh` submits the 500 ps NVT and 5 ns NVE stages for one extracted configuration. `WORK_ROOT` defaults to `./runs`.
