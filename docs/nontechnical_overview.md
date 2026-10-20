# Non-technical project overview

## The project in 30 seconds

Detailed fluid simulations can be extremely expensive because they calculate
the motion of a fluid at millions of locations over thousands of time steps. In
this project, I built a data-driven workflow that learns a much smaller
mathematical model from those simulations. This reduced model is designed to
reproduce the important behaviour of the flow at a fraction of the online
computational cost.

I tested the approach on two difficult examples: a pair of vortices colliding
with a wall, and a heated fluid whose large circulation pattern reverses
direction.

## Why this is useful

High-fidelity simulations are valuable when maximum detail is required, but
they can be too slow when the same system must be evaluated many times, as in
design exploration, optimization, uncertainty analysis, control, or digital
twins.

A reduced-order model acts as a scientific surrogate. It uses data from a
trusted simulator to retain the dominant behaviour in a much smaller system
that can be evaluated more quickly.

## The basic idea

Imagine a flow simulation as a long video in which every frame contains a very
large number of velocity and temperature values. Proper Orthogonal
Decomposition finds recurring spatial patterns in that video, allowing each
frame to be approximated with only a small number of variables.

Operator Inference then learns how those variables evolve over time. Unlike an
unrestricted black-box predictor, it uses a mathematical structure derived from
the governing fluid equations. I also implemented an exact variant that asks
the full simulator a carefully designed set of one-step questions, making it
possible to verify the learned model against a classical physics-based method.

## My contribution

The full-order Navier-Stokes solver was an existing open-source project. I built
the reduced-modelling and experimental layer around it, including:

- the data-reduction and Operator Inference algorithms;
- an independent projection-based implementation for verification;
- global, space-local, and regime-specific reduced models;
- a switched model for the convection-reversal problem;
- automated simulation, analysis, plotting, and experiment-tracking workflows;
- scalable snapshot storage and CPU/GPU execution support;
- physical diagnostics and numerical tests for validation.

The two case studies required different solutions. The vortex-wall collision
creates small, localized structures that challenge a single global compressed
representation. The convection problem changes its overall circulation
direction, so I connected separate models for the periods before, during, and
after the reversal.

## Hard skills developed

- **Julia and scientific software engineering:** modular numerical code,
  configuration-driven experiments, automated tests, and reproducible runs.
- **Data science and scientific machine learning:** dimensionality reduction,
  feature construction, least-squares regression, model selection, and error
  analysis.
- **Numerical linear algebra:** singular-value decomposition, weighted inner
  products, sparse bases, conditioning, and rank-aware algorithms.
- **Computational fluid dynamics:** incompressible Navier-Stokes simulations,
  finite-volume discretization, boundary conditions, and flow diagnostics.
- **High-performance data processing:** large snapshot matrices, memory-mapped
  storage, CPU/GPU backends, and performance benchmarking.
- **Scientific validation:** comparison against independent reference methods,
  physical sanity checks, and transparent analysis of model behaviour.

## Further reading

For the complete methodology, equations, experiments, results, and discussion,
read the [Bachelor thesis (PDF)](../report/Francisco_Abajo_Lis_Bachelor_Thesis.pdf).
The [technical README](../README.md) and
[methodology notes](methodology.md) provide a shorter technical reference.
