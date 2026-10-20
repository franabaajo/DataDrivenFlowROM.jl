# Non-technical project overview

## The project in 30 seconds

Detailed fluid simulations can be extremely expensive because they calculate
the motion of a fluid at millions of locations over thousands of time steps. In
this project, I built a data-driven workflow that learns a much smaller
mathematical model from those simulations. The smaller model aims to reproduce
the important flow behaviour at a fraction of the online computational cost.

I tested the approach on two difficult examples: a pair of vortices colliding
with a wall, and a heated fluid whose large circulation pattern reverses
direction. The work shows both the potential of these reduced models and where
they lose accuracy.

## Why this is useful

A high-fidelity simulation is valuable when maximum detail is required, but it
can be too slow for tasks that require the same system to be evaluated many
times. Examples include design exploration, uncertainty analysis, optimization,
control, and real-time digital twins.

A reduced-order model (ROM) acts as a scientific surrogate. It does not replace
the underlying physics. Instead, it uses data from a trusted simulator to retain
the dominant behaviour in a much smaller system that is cheaper to advance in
time.

## The basic idea

Imagine a flow simulation as a long video in which every frame contains a very
large number of velocity and temperature values. Proper Orthogonal
Decomposition (POD) finds recurring spatial patterns in that video. Each full
frame can then be approximated by combining a relatively small number of these
patterns.

Reducing the number of variables is only half of the problem: the model must
also know how their strengths change with time. Operator Inference learns those
evolution rules from simulation data while retaining a mathematical structure
derived from the governing fluid equations. The result is a compact,
interpretable dynamical system rather than an unrestricted black-box predictor.

I also evaluated an exact form of Operator Inference. Instead of learning only
from one recorded trajectory, it asks the full simulator a carefully designed
set of one-step questions. These samples contain exactly the information needed
to identify the reduced dynamics. This makes it possible to check whether the
learned model matches the classical physics-based projection method to numerical
precision.

## What I built

The full-order Navier-Stokes solver was an existing open-source project. My work
was the reduced-modelling and experimental layer around it, including:

- weighted and centered POD methods adapted to the solver's numerical grid;
- trajectory-based and exact Operator Inference for quadratic flow dynamics;
- an independent intrusive projection implementation for verification;
- global, space-local, and regime-local reduced bases;
- a switched model for the convection-reversal problem;
- automated full-order and reduced-order experiment pipelines;
- streamed, memory-mapped storage for simulation snapshots;
- physical diagnostics for energy, vorticity, heat transport, and circulation;
- reproducible configuration, run tracking, plotting, testing, and CPU/GPU
  execution support.

This required combining numerical linear algebra, regression, fluid mechanics,
scientific software design, data processing, and careful validation.

## The two case studies

### Vortices colliding with a wall

The first case follows a vortex dipole as it travels toward and collides with a
no-slip wall. It is difficult for a reduced model because the collision creates
thin, localized structures near the wall.

The exact data-driven method recovered the same reduced operators as the
intrusive reference to floating-point precision. The reduced model reproduced
the main energetic behaviour, but a single global basis had difficulty
representing the smallest near-wall structures after impact. This distinction
was important: it showed that the dominant remaining error came from the chosen
compressed representation, not from incorrectly learning the dynamics.

### Buoyancy-driven circulation reversal

The second case models a fluid heated from below and cooled from above. The flow
develops a large circulation roll that eventually reverses direction. A single
global reduced model was not sufficient to represent this qualitative change.

I therefore constructed separate local models for the flow before, during, and
after the reversal and connected them with a scheduled switching strategy. The
resulting model captured the reversal, remained stable beyond its training
window, and preserved the sign and approximate strength of the post-reversal
circulation.

## Headline results

- Exact Operator Inference matched intrusive projection to approximately
  `1e-14`-`1e-15` relative operator error.
- Recorded online speedups for selected dipole models were approximately
  `15,100x` and `50,500x`, excluding data generation and model training.
- The selected switched convection model captured the circulation reversal and
  continued stably from the training end at `t = 150` to `t = 300`.
- The experiments identified clear accuracy-cost trade-offs: retaining more
  patterns improves representation but can make quadratic reduced models much
  more expensive.

## What the results mean

The project demonstrates that a physics-structured model can be learned from
simulation queries without assembling the enormous full-order operators. Once
the learned and projected operators agree, the central challenge shifts from
operator identification to choosing a reduced space that contains every flow
regime and localized structure of interest.

It also illustrates an important scientific-computing principle: a fast model
is useful only when its error sources and operating limits are measured. For
that reason, the repository includes negative results and limitations alongside
the strongest accuracy and speed measurements.

## Important limitations

- The reported speedups measure only online time advancement; producing the
  training data and building the ROM have additional costs.
- The circulation switches are scheduled using information from the reference
  simulation, rather than detected autonomously by the reduced model.
- The selected thermal-flow trajectory is a reduced-modelling case study, not a
  grid-converged prediction at the stated physical parameters.
- Very large quadratic reduced models still require techniques such as
  hyper-reduction to remain computationally attractive.

For the equations, detailed numerical results, and reproduction instructions,
return to the [technical README](../README.md) or read the
[methodology and findings](methodology.md).
