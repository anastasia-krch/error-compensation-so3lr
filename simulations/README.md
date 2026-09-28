# Simulation setup

The simulations underlying this work use NPT equilibration, followed by NVT equilibration and NVE production.  The two engines have different input formats, so their complete stage templates are kept separately:

- `ams/` — APPLE&P molecular dynamics with AMS.
- `so3lr/` — SO3LR molecular dynamics with the `so3lr` command-line program.

Each template is parameterised through environment variables. Run NPT once per temperature, select five configurations at 200 ps intervals from its second nanosecond, and use `submit_pipeline.sh` for each 500 ps NVT and 5 ns NVE branch.

The production stage is one uninterrupted NVE run.  Its length can be adjusted with `NSTEPS` (AMS) or `MD_CYCLES` and `MD_STEPS` (SO3LR).

Full periodic configurations at each reported temperature are collected in `representative_frames/`.
