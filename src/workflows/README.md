# Workflow Utilities

Reusable FOM/ROM workflow building blocks for case studies.

## Files

- `Backend.jl`: backend selection and CUDA pressure-solver checks.
- `FOM.jl`: streamed snapshot I/O and FOM stepping helpers.
- `POD.jl`: POD/SVD utility functions.
- `QuadraticROM.jl`: quadratic OpInf fitting and ROM time-integration helpers.

These utilities are case-agnostic and can be included from case-specific scripts under `experiments/`.
