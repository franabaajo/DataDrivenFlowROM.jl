# Rayleigh-Benard convection

This case studies two-dimensional buoyancy-driven convection at
`Ra = 1.6e7`, `Pr = 1`, including a change in direction of the large-scale
circulation.

Centered, energy-weighted POD bases are trained for pre-reversal, transition,
and post-reversal windows. The resulting local affine-quadratic ROMs can be
advanced with scheduled switching and compared across exact OpInf, intrusive
Galerkin projection, and trajectory OpInf.

## Run the demo FOM

From the repository root:

```bash
RAYLEIGH_BACKEND=cpu bash scripts/rayleighbenard.sh run-fom
```

The wrapper defaults to `configs/demo_fom.toml`. Thesis-scale settings are kept
in `configs/fom_base.toml` and require considerably more compute and memory.

## Contents

- `configs/`: demo, full-scale, ROM, and grid-sweep settings;
- `scripts/`: FOM/ROM commands, automatic local-window construction, switching,
  diagnostics, and benchmarks;
- `studies/`: the grid-convergence study.

Reusable algorithms live in `src/workflows/`; thermal diagnostics and
case-specific orchestration live in `src/rayleighbenard/`.
