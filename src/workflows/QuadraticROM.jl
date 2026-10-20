if !isdefined(@__MODULE__, :_WORKFLOW_QUADRATIC_ROM_LOADED)
const _WORKFLOW_QUADRATIC_ROM_LOADED = true

using LinearAlgebra
using ProgressMeter

"""
Return the number of unique quadratic terms `a[j] * a[k]` with `j <= k` for a
state vector of length `r`.
"""
n_unique_quadratic(r::Integer) = div(r * (r + 1), 2)

"""
Return the total number of features for an affine quadratic ROM:

`da/dt = c + L*a + q(a)`.
"""
n_features_affine_quadratic(r::Integer) = 1 + r + n_unique_quadratic(r)

"""
Return the largest ROM dimension `r` such that a quadratic model can be
identified from `K` samples, i.e. `r + r(r+1)/2 <= K`.
"""
function max_quadratic_rank_symmetric(K::Integer)
    r = 0
    while r + div(r * (r + 1), 2) <= K
        r += 1
    end
    return r - 1
end

"""
Return the largest ROM dimension `r` such that an affine quadratic model
(`c + L*a + q(a)`) is identifiable from `K` samples, i.e.
`1 + r + r(r+1)/2 <= K`.
"""
function max_affine_quadratic_rank_symmetric(K::Integer)
    r = 0
    while n_features_affine_quadratic(r + 1) <= K
        r += 1
    end
    return r
end

"""
Return index pairs `(j, k)` for all unique quadratic terms with `1 <= j <= k <= r`.
"""
function quadratic_pairs(r::Integer)
    pairs = Tuple{Int,Int}[]
    for j in 1:r
        for k in j:r
            push!(pairs, (j, k))
        end
    end
    return pairs
end

"""
Build the symmetric quadratic feature matrix for ROM coordinates.

Given `a` with shape `(r, K)` (columns are snapshots), return:
- `Φ` with shape `(K, r + n_unique_quadratic(r))`
- `pairs` describing which `(j, k)` produced each quadratic column.
"""
function quadratic_features_symmetric(a::AbstractMatrix)
    r, K = size(a)
    pairs = quadratic_pairs(r)
    nq = length(pairs)

    T = eltype(a)
    Φ = Matrix{T}(undef, K, r + nq)

    @inbounds for n in 1:K
        an = @view a[:, n]
        Φ[n, 1:r] = an
        for (m, (j, k)) in enumerate(pairs)
            Φ[n, r + m] = an[j] * an[k]
        end
    end

    return Φ, pairs
end

"""
Build the symmetric affine quadratic feature matrix for ROM coordinates.

Given `a` with shape `(r, K)` (columns are snapshots), return:
- `Φ` with shape `(K, 1 + r + n_unique_quadratic(r))` (first column is ones)
- `pairs` describing which `(j, k)` produced each quadratic column.
"""
function quadratic_features_affine_symmetric(a::AbstractMatrix)
    r, K = size(a)
    pairs = quadratic_pairs(r)
    nq = length(pairs)

    T = eltype(a)
    Φ = Matrix{T}(undef, K, 1 + r + nq)

    @inbounds for n in 1:K
        an = @view a[:, n]
        Φ[n, 1] = one(T)
        Φ[n, 2:(1 + r)] = an
        for (m, (j, k)) in enumerate(pairs)
            Φ[n, 1 + r + m] = an[j] * an[k]
        end
    end

    return Φ, pairs
end

"""
DISCLAIMER: No guarantee of stability or accuracy of the fitted ROM. Function not used in Exact OpInf method.
Fit a quadratic ROM `da/dt = L*a + q(a)` from trajectory data.

`a` and `dadt` are `(r, K)`. If `lambda > 0`, regularization is used.
Returns fitted operators and diagnostics as a named tuple.
"""
function fit_quadratic_rom_symmetric(a::AbstractMatrix, dadt::AbstractMatrix; lambda::Real = 0.1)
    r, K = size(a)
    Φ, pairs = quadratic_features_symmetric(a)
    Y = Matrix(dadt')

    Ξ = if lambda == 0
        Φ \ Y
    else
        (Φ' * Φ + lambda * I) \ (Φ' * Y)
    end

    L = Matrix(Ξ[1:r, :]')
    B = Matrix(Ξ[r+1:end, :]')

    Q = zeros(eltype(a), r, r, r)
    @inbounds for i in 1:r
        for (m, (j, k)) in enumerate(pairs)
            c = B[i, m]
            if j == k
                Q[i, j, k] = c
            else
                Q[i, j, k] = c / 2
                Q[i, k, j] = c / 2
            end
        end
    end

    σ = svdvals(Φ)
    return (; L, Q, B, Φ, σ, pairs)
end

"""
DISCLAIMER: No guarantee of stability or accuracy of the fitted ROM.
Fit an affine quadratic ROM `da/dt = c + L*a + q(a)` from trajectory data.
"""
function fit_affine_quadratic_rom_symmetric(a::AbstractMatrix, dadt::AbstractMatrix; lambda::Real = 0.1)
    r, K = size(a)
    Φ, pairs = quadratic_features_affine_symmetric(a)
    Y = Matrix(dadt')

    Ξ = if lambda == 0
        Φ \ Y
    else
        (Φ' * Φ + lambda * I) \ (Φ' * Y)
    end

    c = Vector(Ξ[1, :])
    L = Matrix(Ξ[2:(1 + r), :]')
    B = Matrix(Ξ[(2 + r):end, :]')

    Q = zeros(eltype(a), r, r, r)
    @inbounds for i in 1:r
        for (m, (j, k)) in enumerate(pairs)
            coeff = B[i, m]
            if j == k
                Q[i, j, k] = coeff
            else
                Q[i, j, k] = coeff / 2
                Q[i, k, j] = coeff / 2
            end
        end
    end

    σ = svdvals(Φ)
    return (; c, L, Q, B, Φ, σ, pairs)
end

"""
Construct a deterministic set of ROM states that spans all linear and unique
quadratic features for dimension `r`.
Returns a matrix of shape `(r, r + n_unique_quadratic(r))` where the first `r` columns are the standard basis and the remaining columns are sums of pairs of
basis vectors corresponding to the unique quadratic features. This ensures the feature matrix is full rank for exact OpInf fitting.
"""
function rank_ensuring_states_quadratic(r::Integer; T::Type = Float64)
    nq = n_unique_quadratic(r)
    K = r + nq
    A = zeros(T, r, K)
    col = 1

    for j in 1:r
        A[j, col] = one(T)
        col += 1
    end

    for (j, k) in quadratic_pairs(r)
        A[j, col] = one(T)
        A[k, col] += one(T)
        col += 1
    end

    return A
end

"""
Construct a deterministic set of ROM states that spans constant, linear, and
unique quadratic features for an affine quadratic ROM.

Returns a matrix of shape `(r, 1 + r + n_unique_quadratic(r))` where:
- column 1 is the zero state (captures the constant feature)
- next `r` columns are the standard basis
- remaining columns are sums of pairs of basis vectors corresponding to the
  unique quadratic features.
"""
function rank_ensuring_states_affine_quadratic(r::Integer; T::Type = Float64)
    nq = n_unique_quadratic(r)
    K = 1 + r + nq
    A = zeros(T, r, K)
    col = 2

    for j in 1:r
        A[j, col] = one(T)
        col += 1
    end

    for (j, k) in quadratic_pairs(r)
        A[j, col] = one(T)
        A[k, col] += one(T)
        col += 1
    end

    return A
end

"""
Estimate the exact-OpInf one-step training size following the MATLAB reference
formula used in `dt_estimate.m` for the no-input quadratic case (`is = [1, 2]`).

Given reduced snapshots `a` and reduced derivatives `dadt`, both `(r, K)`, this
uses only the first mode:

`dt = 1 / max_k |dadt₁(k)| / ‖[a₁(k), a₁(k)^2]‖₂`.

As in the MATLAB script, only columns `1:(K-1)` are used when `K > 1`.

If the ratio is non-finite or zero for all samples, falls back to `fallback_dt`.
"""
function estimate_exact_opinf_dt_matlab(
    a::AbstractMatrix,
    dadt::AbstractMatrix;
    fallback_dt::Real,
)
    size(a) == size(dadt) || error("Size mismatch: a has size $(size(a)) but dadt has size $(size(dadt))")
    size(a, 1) >= 1 || error("Need at least one reduced mode to estimate dt")
    size(a, 2) >= 1 || error("Need at least one snapshot to estimate dt")

    kend = size(a, 2) > 1 ? (size(a, 2) - 1) : 1
    a1 = @view a[1, 1:kend]
    d1 = @view dadt[1, 1:kend]

    maxratio = 0.0
    @inbounds for k in eachindex(a1)
        ak = Float64(a1[k])
        d1k = Float64(d1[k])
        denom = hypot(ak, ak * ak) # Equivalent to vecwise_2norm([a1; a1^2]).
        ratio = if denom == 0
            d1k == 0 ? 0.0 : Inf
        else
            abs(d1k) / denom
        end
        if isfinite(ratio)
            maxratio = max(maxratio, ratio)
        else
            maxratio = Inf
            break
        end
    end

    if !isfinite(maxratio) || maxratio <= 0
        fallback = Float64(fallback_dt)
        @warn "MATLAB-style dt estimate was non-finite or zero; using fallback dt." maxratio fallback
        return fallback
    end
    return 1.0 / maxratio
end

"""
Fit a quadratic ROM in the exact OpInf setting from feature-complete training
states `Atrain` and derivatives `Dtrain`.

Returns fitted operators and residual metrics (`resnorm`, `relres`).
"""
function fit_quadratic_rom_exact_symmetric(Atrain::AbstractMatrix, Dtrain::AbstractMatrix)
    r, K = size(Atrain)
    Φ, pairs = quadratic_features_symmetric(Atrain)
    Y = Matrix(Dtrain')

    K == size(Φ, 1) || error("Feature matrix size mismatch: K = $K, Φ rows = $(size(Φ, 1))")
    size(Φ, 2) == r + n_unique_quadratic(r) ||
        error("Feature matrix has wrong column count for r = $r")
    if size(Φ, 1) != size(Φ, 2)
        @warn "Exact OpInf feature matrix is not square; check training state construction." sizeΦ = size(Φ)
    end

    Ξ = Φ \ Y # No regularization since features are designed to be full rank.

    L = Matrix(Ξ[1:r, :]')
    B = Matrix(Ξ[r+1:end, :]') # Coefficients for unique quadratic terms; need to symmetrize into Q tensor.

    Q = zeros(eltype(Atrain), r, r, r)
    @inbounds for i in 1:r
        for (m, (j, k)) in enumerate(pairs)
            c = B[i, m]
            if j == k
                Q[i, j, k] = c
            else
                Q[i, j, k] = c / 2
                Q[i, k, j] = c / 2
            end
        end
    end

    σ = svdvals(Φ)
    residual = Φ * Ξ - Y # Should be zero if features are full rank and no regularization. Calculated for diagnostics (machine precision).
    resnorm = norm(residual) # Absolute residual norm.
    relres = resnorm / max(norm(Y), eps(eltype(residual))) # Relative residual norm.
    return (; L, Q, B, Φ, σ, pairs, resnorm, relres)
end

"""
Fit an affine quadratic ROM in the exact OpInf setting:

`da/dt = c + L*a + q(a)`.
"""
function fit_affine_quadratic_rom_exact_symmetric(Atrain::AbstractMatrix, Dtrain::AbstractMatrix)
    r, K = size(Atrain)
    Φ, pairs = quadratic_features_affine_symmetric(Atrain)
    Y = Matrix(Dtrain')

    K == size(Φ, 1) || error("Feature matrix size mismatch: K = $K, Φ rows = $(size(Φ, 1))")
    size(Φ, 2) == n_features_affine_quadratic(r) ||
        error("Feature matrix has wrong column count for r = $r")
    if size(Φ, 1) != size(Φ, 2)
        @warn "Exact OpInf affine feature matrix is not square; check training state construction." sizeΦ = size(Φ)
    end

    Ξ = Φ \ Y

    c = Vector(Ξ[1, :])
    L = Matrix(Ξ[2:(1 + r), :]')
    B = Matrix(Ξ[(2 + r):end, :]')

    Q = zeros(eltype(Atrain), r, r, r)
    @inbounds for i in 1:r
        for (m, (j, k)) in enumerate(pairs)
            coeff = B[i, m]
            if j == k
                Q[i, j, k] = coeff
            else
                Q[i, j, k] = coeff / 2
                Q[i, k, j] = coeff / 2
            end
        end
    end

    σ = svdvals(Φ)
    residual = Φ * Ξ - Y
    resnorm = norm(residual)
    relres = resnorm / max(norm(Y), eps(eltype(residual)))
    return (; c, L, Q, B, Φ, σ, pairs, resnorm, relres)
end

"""
Evaluate the quadratic ROM right-hand side in place.

Overwrites `da` with `L*a + q(a)` and returns `da`.
"""
function rom_rhs!(da, a, L, Q)
    mul!(da, L, a)
    r = length(a)
    @inbounds for i in 1:r
        s = zero(eltype(a))
        for k in 1:r, j in 1:r
            s += Q[i, j, k] * a[j] * a[k]
        end
        da[i] += s
    end
    return da
end

"""
Evaluate the affine quadratic ROM right-hand side in place.

Overwrites `da` with `c + L*a + q(a)` and returns `da`.
"""
function rom_rhs!(da, a, c, L, Q)
    copyto!(da, c)
    T = eltype(a)
    mul!(da, L, a, one(T), one(T))
    r = length(a)
    @inbounds for i in 1:r
        s = zero(T)
        for k in 1:r, j in 1:r
            s += Q[i, j, k] * a[j] * a[k]
        end
        da[i] += s
    end
    return da
end

"""
Evaluate the Jacobian of the quadratic ROM in place at state `a`.

Overwrites `J` with `∂(L*a + q(a))/∂a` and returns `J`.
"""
function rom_jacobian!(J, a, L, Q)
    copyto!(J, L)
    r = length(a)
    @inbounds for i in 1:r, m in 1:r
        s = zero(eltype(a))
        for k in 1:r
            s += Q[i, m, k] * a[k] + Q[i, k, m] * a[k]
        end
        J[i, m] += s
    end
    return J
end

"""
Propose a stable ROM time step from the Jacobian 1-norm.

Returns `min(dtmax, safety / max(opnorm(J, 1), 1e-12))`.
"""
function propose_rom_timestep(a, L, Q; dtmax::Real, safety::Real = 0.15)
    r = length(a)
    J = Matrix{eltype(a)}(undef, r, r)
    rom_jacobian!(J, a, L, Q)
    ρ = opnorm(J, 1)
    return min(dtmax, safety / max(ρ, 1e-12))
end

"""
Advance the ROM state by one RK4 step of size `dt` in place.

Uses preallocated work arrays `k1`, `k2`, `k3`, `k4`, and `atmp`.
"""
function rk4_step!(a, dt, L, Q, k1, k2, k3, k4, atmp)
    rom_rhs!(k1, a, L, Q)

    @. atmp = a + 0.5 * dt * k1
    rom_rhs!(k2, atmp, L, Q)

    @. atmp = a + 0.5 * dt * k2
    rom_rhs!(k3, atmp, L, Q)

    @. atmp = a + dt * k3
    rom_rhs!(k4, atmp, L, Q)

    @. a = a + dt * (k1 + 2k2 + 2k3 + k4) / 6
    return a
end

"""
Advance the affine quadratic ROM state by one RK4 step of size `dt` in place.
"""
function rk4_step!(a, dt, c, L, Q, k1, k2, k3, k4, atmp)
    rom_rhs!(k1, a, c, L, Q)

    @. atmp = a + 0.5 * dt * k1
    rom_rhs!(k2, atmp, c, L, Q)

    @. atmp = a + 0.5 * dt * k2
    rom_rhs!(k3, atmp, c, L, Q)

    @. atmp = a + dt * k3
    rom_rhs!(k4, atmp, c, L, Q)

    @. a = a + dt * (k1 + 2k2 + 2k3 + k4) / 6
    return a
end

"""
Simulate a quadratic ROM with adaptive substepping and snapshot output.

Integrates from `a0` to `tsim`, stores `nsave + 1` states (including the
initial condition), and rejects unstable trial steps by halving `dt`.
"""
function simulate_rom_adaptive(
    a0,
    L,
    Q;
    tsim::Real,
    nsave::Integer,
    safety::Real = 0.03,
    growthfac::Real = 1.05,
    dtmin::Real = 1e-10,
)
    r = length(a0)
    T = eltype(a0)

    a = copy(a0)
    atrial = similar(a)
    A = zeros(T, r, nsave + 1)
    copyto!(view(A, :, 1), a)

    k1 = similar(a)
    k2 = similar(a)
    k3 = similar(a)
    k4 = similar(a)
    atmp = similar(a)
    Jtmp = Matrix{T}(undef, r, r)

    t = zero(T)
    prog = Progress(nsave; desc = "Running ROM")

    for isnap in 1:nsave
        tstop = isnap / nsave * tsim
        isubstep = 0

        while t < prevfloat(tstop)
            dtcap = tstop - t
            dt = min(propose_rom_timestep(a, L, Q; dtmax = dtcap, safety = safety), dtcap)

            accepted = false
            while !accepted
                copyto!(atrial, a)
                rk4_step!(atrial, dt, L, Q, k1, k2, k3, k4, atmp)

                bad = !all(isfinite, atrial)
                na = norm(a)
                ntrial = norm(atrial)
                growth_bad = ntrial > growthfac * max(na, 1e-12)

                if bad || growth_bad
                    rom_jacobian!(Jtmp, a, L, Q)
                    @info "Rejected ROM step" (
                        t = t,
                        dt = dt,
                        jacnorm = opnorm(Jtmp, 1),
                        norma = na,
                        normtrial = ntrial,
                        bad = bad,
                        growth_bad = growth_bad,
                    )

                    dt *= 0.5
                    if dt < dtmin
                        error("ROM timestep became too small at t = $t")
                    end
                else
                    copyto!(a, atrial)
                    t += dt
                    isubstep += 1
                    accepted = true
                end
            end
        end

        copyto!(view(A, :, isnap + 1), a)
        next!(prog; showvalues = [("substeps", isubstep), ("time", t)])
    end

    return A
end

"""
Simulate an affine quadratic ROM `da/dt = c + L*a + q(a)` with adaptive substepping.
"""
function simulate_rom_adaptive(
    a0,
    c,
    L,
    Q;
    tsim::Real,
    nsave::Integer,
    safety::Real = 0.03,
    growthfac::Real = 1.05,
    dtmin::Real = 1e-10,
)
    r = length(a0)
    T = eltype(a0)
    length(c) == r || error("Constant term length mismatch: length(c) = $(length(c)) but length(a0) = $r")

    a = copy(a0)
    atrial = similar(a)
    A = zeros(T, r, nsave + 1)
    copyto!(view(A, :, 1), a)

    k1 = similar(a)
    k2 = similar(a)
    k3 = similar(a)
    k4 = similar(a)
    atmp = similar(a)
    Jtmp = Matrix{T}(undef, r, r)

    t = zero(T)
    prog = Progress(nsave; desc = "Running ROM")

    for isnap in 1:nsave
        tstop = isnap / nsave * tsim
        isubstep = 0

        while t < prevfloat(tstop)
            dtcap = tstop - t
            dt = min(propose_rom_timestep(a, L, Q; dtmax = dtcap, safety = safety), dtcap)

            accepted = false
            while !accepted
                copyto!(atrial, a)
                rk4_step!(atrial, dt, c, L, Q, k1, k2, k3, k4, atmp)

                bad = !all(isfinite, atrial)
                na = norm(a)
                ntrial = norm(atrial)
                growth_bad = ntrial > growthfac * max(na, 1e-12)

                if bad || growth_bad
                    rom_jacobian!(Jtmp, a, L, Q)
                    @info "Rejected ROM step" (
                        t = t,
                        dt = dt,
                        jacnorm = opnorm(Jtmp, 1),
                        norma = na,
                        normtrial = ntrial,
                        bad = bad,
                        growth_bad = growth_bad,
                    )

                    dt *= 0.5
                    if dt < dtmin
                        error("ROM timestep became too small at t = $t")
                    end
                else
                    copyto!(a, atrial)
                    t += dt
                    isubstep += 1
                    accepted = true
                end
            end
        end

        copyto!(view(A, :, isnap + 1), a)
        next!(prog; showvalues = [("substeps", isubstep), ("time", t)])
    end

    return A
end

end
