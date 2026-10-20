# Reproducing the experiments

## Requirements

- Julia 1.10 or newer;
- Git;
- a POSIX-compatible shell for the convenience wrappers, or PowerShell with
  the direct Julia commands;
- optional NVIDIA GPU and compatible CUDA driver for full-scale runs.

Install dependencies from the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Lightweight verification

Run the tests and the small CPU dipole example:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
DIPOLE_BACKEND=cpu bash scripts/dipole.sh run-all
```

The demo verifies the pipeline, not the thesis error values. Its coarse grid
cannot resolve the wall collision used in the scientific comparison.

## Thesis-scale dipole study

The full configuration is
`experiments/DipoleWallCollision/configs/fom_base.toml`. It uses a `1536 x 1536`
grid and is computationally expensive.

```bash
DIPOLE_BACKEND=gpu bash scripts/dipole.sh run-fom \
  experiments/DipoleWallCollision/configs/fom_base.toml
```

Use the printed FOM run ID to execute a ROM:

```bash
bash scripts/dipole.sh run-rom FOM_RUN_ID \
  experiments/DipoleWallCollision/configs/rom_base.toml
```

## Rayleigh-Benard study

Generate the FOM and record its run ID:

```bash
RAYLEIGH_BACKEND=gpu bash scripts/rayleighbenard.sh run-fom \
  experiments/RayleighBenardConvection/configs/fom_base.toml
```

Then run a single ROM using:

```bash
bash scripts/rayleighbenard.sh run-rom FOM_RUN_ID \
  experiments/RayleighBenardConvection/configs/rom_base.toml
```

The advanced switched-ROM commands are documented by:

```bash
julia --project=. experiments/RayleighBenardConvection/scripts/rayleigh_tasks.jl --help
```

## Outputs

Generated data are placed under `output/<case>/` and include:

- exact TOML configuration;
- metadata and timing information;
- streamed snapshot matrices;
- fitted basis and reduced operators;
- reduced trajectories;
- figures and diagnostic tables;
- a registry row linking each ROM to its parent FOM.

## Reproducing the thesis tables

The scripts under `analysis/` expect the original full-scale run directories.
Because those binary datasets are too large for Git, the final compact tables
are retained under `results/tables/`. Update the run IDs at the top of an
analysis script if regenerated runs receive different identifiers.

