# Space-local POD

## Motivation

A global POD basis minimizes snapshot projection error at a fixed rank, but its
basis vectors are dense over the entire domain. For large flow states this makes
projection and reconstruction expensive and can poorly emphasize localized
structures such as boundary layers.

Space-local POD partitions the state degrees of freedom into non-overlapping
spatial blocks and computes a POD basis independently in each block. The
assembled basis is block sparse.

## Formulation

For a snapshot matrix `X` with row blocks `I_1, ..., I_C`, define

```text
X_i = X[I_i, :].
```

After computing a weighted local SVD for each `X_i`, retain `r_i` modes and
assemble the global basis from the local modes. The total rank is

```text
r_total = sum(r_i).
```

Because the blocks do not overlap, the global projection error decomposes into
the sum of discarded local singular-value energies. With finite-volume row
weights `w`, the implementation scales local rows by `sqrt(w)` and constructs a
basis satisfying

```text
Phi' * Diagonal(w) * Phi approximately I.
```

## Rank-selection strategies

- `S1`: pool local singular values and retain the globally largest values while
  keeping a fixed total rank.
- `S3`: choose each block's rank independently using a local energy threshold
  or explicit local ranks.

Partitions can be uniform in `x`, uniform in `x-y`, or defined by manual cut
locations. Velocity and scalar degrees of freedom are assigned using their
staggered-grid coordinates.

## Trade-offs

At the same total rank, global POD has the lowest possible projection error.
Space-local POD is useful because it produces sparse, localized modes and can
reduce the cost of compatible online operations. It may require a larger total
rank for comparable projection accuracy.

For quadratic Operator Inference, rank inflation is especially important:
the number of unique quadratic features grows as `r(r+1)/2`. A partition that
improves spatial locality can therefore worsen identifiability or dense ROM
cost unless the nonlinear operators also exploit sparsity or hyper-reduction.

The implementation is in `src/workflows/POD.jl` and its weighted
orthonormality and block support are checked in
`test/workflow_pod_space_local.jl`.

