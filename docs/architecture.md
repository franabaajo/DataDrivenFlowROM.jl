# Software architecture

## Reusable layer

`src/workflows/` contains case-independent building blocks:

- `FOM.jl`: streamed snapshots, memory mapping, and one-step FOM queries;
- `POD.jl`: weighted Gram-matrix POD and space-local sparse bases;
- `QuadraticROM.jl`: feature construction, operator fitting, rank-ensuring
  states, adaptive RK4 integration, and stability checks;
- `Backend.jl`: CPU/CUDA selection used by the experiment runners.

The package entry point exposes the reusable numerical algorithms. Backend and
plotting dependencies remain in the experiment layer so the mathematical core
stays understandable in isolation.

## Case layer

Each directory under `src/dipole/` and `src/rayleighbenard/` separates:

- typed TOML configuration;
- geometry and initial conditions;
- deterministic paths and run identifiers;
- FOM execution and snapshot loading;
- ROM construction and simulation;
- diagnostic calculations;
- plotting and registry management.

The Julia entry points under `experiments/*/scripts/` compose those files into
command-line workflows. This keeps reusable algorithms independent of
case-specific physics while preserving transparent research scripts.

## Data lifecycle

```text
TOML configuration
      -> deterministic run ID
      -> streamed FOM snapshots + metadata
      -> basis and reduced operators
      -> ROM trajectory
      -> diagnostics, plots, and registry entry
```

Large binary snapshots are reproducible outputs and are ignored by Git. Compact
CSV summaries and selected publication figures are committed.

