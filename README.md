# DataDrivenFlowROM.jl

[![CI](https://github.com/franabaajo/DataDrivenFlowROM.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/franabaajo/DataDrivenFlowROM.jl/actions/workflows/ci.yml)
[![Julia 1.10](https://img.shields.io/badge/Julia-1.10%2B-9558B2.svg)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Data-driven reduced-order modelling of two-dimensional incompressible turbulent
flows in Julia.

This project implements and evaluates Proper Orthogonal Decomposition (POD),
trajectory-based Operator Inference, exact Operator Inference, and intrusive
Galerkin projection on two challenging flow problems:

- collision of a vortex dipole with a no-slip wall at `Re = 2500`;
- Rayleigh-Benard convection at `Ra = 1.6e7`, including a reversal of the
  large-scale circulation.

It was developed for my combined Applied Mathematics and Applied Physics
Bachelor End Project at Eindhoven University of Technology. Full-order
simulations use the MIT-licensed
[IncompressibleNavierStokes.jl](https://github.com/agdestein/IncompressibleNavierStokes.jl)
solver. The ROM algorithms, experiment automation, diagnostics, analysis, and
case studies in this repository are my contribution.

![Rayleigh-Benard circulation reversal](results/figures/rayleigh_fom_reversal_snapshots.png)

## Highlights

- Energy-weighted POD consistent with the solver's staggered finite-volume grid.
- Centered POD and fluctuation-energy rank selection for thermal flows.
- Affine-quadratic and quadratic reduced dynamical systems.
- Exact OpInf through rank-ensuring one-step full-order queries.
- Matrix-free intrusive Galerkin operators for verification.
- Space-local POD with sparse bases and configurable domain partitions.
- Scheduled local ROMs for pre-reversal, transition, and post-reversal regimes.
- Streamed, memory-mapped snapshots for large simulations.
- CPU/GPU selection, TOML configurations, deterministic run identifiers, and
  experiment registries.
- Reproducible physical diagnostics: kinetic energy, enstrophy, Nusselt number,
  and normalized angular momentum.

## Key results

| Experiment | Result |
|---|---|
| Dipole exact OpInf | Intrusive linear and quadratic operators recovered with relative errors of `2.30e-13` and `1.77e-14` |
| Rayleigh-Benard exact OpInf | Local affine-quadratic operators recovered to approximately `1e-14`-`1e-15` relative error |
| Dipole online performance | Recorded median speedups of about `15,100x` at rank 22 and `50,500x` at rank 13 |
| Circulation reversal | The `eta = 0.99` switched model achieved mean relative state error `0.16034` and angular-momentum RMSE `0.07033` |
| Post-training prediction | Stable continuation from training end `t = 150` to `t = 300`, with mean post-training state error `0.18786` |

The timings above measure online advancement and exclude full-order data
generation. See [Scientific findings](docs/methodology.md) for interpretation
and limitations.

![Exact and intrusive operators produce matching switched-ROM errors](results/figures/rayleigh_method_relative_error.png)

## Repository map

```text
src/workflows/                 reusable POD, FOM-query, and ROM algorithms
src/dipole/                    dipole configuration, simulation, ROM, and diagnostics
src/rayleighbenard/            thermal-flow configuration, ROM, switching, and diagnostics
experiments/                   executable case studies and TOML configurations
analysis/                      scripts used to produce the reported statistics
results/                       curated figures and compact result tables
test/                          contribution-focused numerical and I/O tests
report/                        bachelor thesis
docs/                          methodology and full reproducibility instructions
```

## Quick start

Install Julia 1.10 or newer, clone the repository, and instantiate its
environment:

```bash
git clone https://github.com/franabaajo/DataDrivenFlowROM.jl.git
cd DataDrivenFlowROM.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run the lightweight CPU dipole pipeline:

```bash
DIPOLE_BACKEND=cpu bash scripts/dipole.sh run-all
```

On Windows PowerShell, the equivalent commands are:

```powershell
$env:DIPOLE_BACKEND = "cpu"
julia --project=. experiments/DipoleWallCollision/scripts/dipole_tasks.jl run-all `
  experiments/DipoleWallCollision/configs/demo_fom.toml `
  experiments/DipoleWallCollision/configs/demo_rom.toml
```

The demo uses a deliberately small grid as a smoke test. Thesis-scale
configurations are retained as `fom_base.toml` and `rom_base.toml`; they require
substantially more memory and compute time.

For Rayleigh-Benard convection, first generate a FOM run:

```bash
RAYLEIGH_BACKEND=cpu bash scripts/rayleighbenard.sh run-fom
```

The command prints the deterministic run identifier needed by `run-rom`.

## Methods at a glance

```text
Full-order simulation
        |
        v
weighted snapshots --> POD / centered POD --> reduced coordinates
                                                |
                    +---------------------------+------------------+
                    |                           |                  |
                    v                           v                  v
          trajectory OpInf             exact OpInf       intrusive Galerkin
                    |                           |                  |
                    +---------------------------+------------------+
                                                |
                                                v
                              adaptive reduced-order simulation
                                                |
                                                v
                          physical diagnostics and FOM comparison
```

Exact OpInf and intrusive projection agree when they use the same basis and
polynomial feature space. That agreement does not imply exact recovery of the
full-order trajectory: POD truncation, regime coverage, and local structures
remain important sources of error.

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The tests cover weighted space-local POD, snapshot-file integrity, conductive
Nusselt-number sanity checks, and Rayleigh-Benard training-window behavior.

## Report and reproducibility

- [Bachelor thesis](report/Francisco_Abajo_Lis_Bachelor_Thesis.pdf)
- [Methodology and findings](docs/methodology.md)
- [Space-local POD](docs/space_local_pod.md)
- [Reproducing the experiments](docs/reproducibility.md)
- [Software architecture](docs/architecture.md)

Large full-order snapshots are intentionally not versioned. The repository
contains the configurations and scripts needed to regenerate them, while the
compact numerical summaries used in the thesis are retained under `results/`.

## Scope and limitations

- The `256 x 256` Rayleigh-Benard trajectory used for the switched-ROM study is
  a selected reduced-modelling trajectory, not a grid-converged prediction at
  the stated Rayleigh number.
- Switching times are derived from the training FOM; the demonstrated model is
  scheduled rather than fully autonomous.
- High-rank quadratic ROMs remain expensive without hyper-reduction.
- Exact operator recovery removes identification error relative to intrusive
  projection, but not reduced-basis truncation or closure error.

## Authorship and attribution

Copyright (c) 2026 Francisco Abajo Lis.

This repository is licensed under the MIT License. The full-order solver is a
separate upstream project with its own authors and license. See
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) and `LICENSE-UPSTREAM`.

