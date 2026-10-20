using Dates
using JLD2
using KernelAbstractions
using LinearAlgebra
import IncompressibleNavierStokes as NS

include(joinpath(@__DIR__, "..", "workflows", "QuadraticROM.jl"))

function _pod_energy_retained(S, r::Integer)
    isempty(S) && return 0.0
    r < 1 && return 0.0
    r_eff = min(Int(r), length(S))
    tot = sum(abs2, S)
    tot <= 0 && return 0.0
    return sum(abs2, @view(S[1:r_eff])) / tot
end

function centered_pod_spectrum_from_weighted_gram(
    X::AbstractMatrix,
    w::AbstractVector;
    chunk_rows::Integer = 25_000,
)
    N, K = size(X)
    length(w) == N || error("Weight length mismatch: length(w) = $(length(w)) but size(X, 1) = $N")
    chunk_rows >= 1 || error("chunk_rows must be >= 1, got $chunk_rows")

    T = promote_type(eltype(X), eltype(w))
    state_mean = vec(sum(X; dims = 2)) ./ K
    gram = zeros(T, K, K)

    row = 1
    while row <= N
        lastrow = min(row + Int(chunk_rows) - 1, N)
        Xc = Matrix{T}(@view X[row:lastrow, :])
        Xc .-= reshape(@view(state_mean[row:lastrow]), :, 1)
        Y = sqrt.(T.(@view w[row:lastrow])) .* Xc
        mul!(gram, transpose(Y), Y, one(T), one(T))
        row = lastrow + 1
    end

    eig = eigen(Symmetric(gram))
    order = sortperm(eig.values; rev = true)
    values = max.(eig.values[order], zero(eltype(eig.values)))
    S = sqrt.(values)
    V = eig.vectors[:, order]
    return (; S, V, state_mean)
end

function centered_pod_mode_from_gram(X::AbstractMatrix, state_mean::AbstractVector, spectrum, imode::Integer)
    1 <= imode <= length(spectrum.S) || error("POD mode index out of bounds: $imode")
    sigma = spectrum.S[imode]
    sigma_max = isempty(spectrum.S) ? zero(sigma) : maximum(spectrum.S)
    tol = sqrt(eps(Float64)) * max(Float64(sigma_max), 1.0)
    if Float64(sigma) <= tol
        return zeros(eltype(X), size(X, 1))
    end
    v = @view spectrum.V[:, imode]
    mode = X * v
    mode .-= state_mean .* sum(v)
    mode ./= sigma
    return mode
end

function project_snapshots_with_mean(projector::AbstractMatrix, X::AbstractMatrix, state_mean)
    A = projector * X
    if !isnothing(state_mean)
        amean = projector * state_mean
        A .-= reshape(amean, :, 1)
    end
    return A
end

function full_state_from_reduced(rom_basis::AbstractMatrix, a::AbstractVector, state_mean)
    x = rom_basis * a
    isnothing(state_mean) || (x .+= state_mean)
    return x
end

"""
Compute the ROM basis and projector from FOM snapshots.

Returns a named tuple with `rom_basis`, `projector`, singular values, and
selected basis indices.
"""
function build_rom_basis(
    snapshots_state::AbstractMatrix,
    setup,
    romcfg::ROMConfig;
    requested_rank::Int,
    snapshot_indices = 1:size(snapshots_state, 2),
)
    Kfull = size(snapshots_state, 2)
    length(snapshot_indices) == Kfull ||
        error("Basis snapshot-index count mismatch: got $(length(snapshot_indices)) indices for $Kfull snapshots")
    idx_basis = snapshot_indices
    Xbasis = snapshots_state

    # Use a weighted inner product that matches the discrete energy norm on
    # stretched/staggered grids. This avoids over-emphasizing regions with finer
    # mesh spacing (e.g. near walls).
    w = state_mass_vector(setup, NS.vectorfield(setup), NS.scalarfield(setup))
    length(w) == size(Xbasis, 1) || error("State mass vector length mismatch with snapshot dimension")

    basis_mode = lowercase(romcfg.basis_mode)
    if basis_mode == "global"
        rank_sel_raw = lowercase(strip(romcfg.rank_selection))
        rank_selection = if rank_sel_raw in ("npod", "fixed", "rank")
            "npod"
        elseif rank_sel_raw in ("energy_fraction", "energyfraction", "energy", "fraction", "eta")
            "energy_fraction"
        else
            error("Unsupported rank_selection: $(romcfg.rank_selection). Use \"npod\" or \"energy_fraction\".")
        end
        use_centered_basis = romcfg.basis_centered
        raw_spectrum = pod_spectrum_from_weighted_gram(Xbasis, w)
        centered_spectrum = use_centered_basis ?
            centered_pod_spectrum_from_weighted_gram(Xbasis, w) : nothing
        spectrum = use_centered_basis ? centered_spectrum : raw_spectrum

        rfit = if rank_selection == "npod"
            r = min(requested_rank, Kfull)
            if r < requested_rank
                @warn "Requested npod = $(requested_rank), but only $Kfull snapshots are available. Using rfit = $r."
            end
            r
        else
            ef = romcfg.energy_fraction
            (0.0 < ef <= 1.0) || error("energy_fraction must be in (0, 1] for rank_selection=\"energy_fraction\", got $ef")
            r = rank_from_energy_safe(spectrum.S; energy_fraction = ef)
            r >= 1 || error("Energy-based rank selection returned r = $r; check snapshot energy content.")
            r
        end

        energy_retained = _pod_energy_retained(spectrum.S, rfit)
        raw_energy_retained = _pod_energy_retained(raw_spectrum.S, rfit)
        fluctuation_energy_retained = use_centered_basis ?
            energy_retained : nothing
        state_mean = use_centered_basis ? spectrum.state_mean : nothing
        rom_basis = Matrix{eltype(Xbasis)}(undef, size(Xbasis, 1), rfit)
        for i in 1:rfit
            mode = use_centered_basis ?
                centered_pod_mode_from_gram(Xbasis, state_mean, spectrum, i) :
                pod_mode_from_gram(Xbasis, spectrum, i)
            copyto!(@view(rom_basis[:, i]), mode)
        end
        projector = rom_basis' * Diagonal(w)
        rank_info = (;
            rank_selection = rank_selection,
            basis_centered = use_centered_basis,
            energy_basis = use_centered_basis ? "centered_fluctuation" : "raw_state",
            requested_npod = requested_rank,
            energy_fraction = romcfg.energy_fraction,
            energy_retained = energy_retained,
            raw_energy_retained = raw_energy_retained,
            fluctuation_energy_retained = fluctuation_energy_retained,
            rfit = rfit,
        )
        return (; rom_basis, projector, singular_values = spectrum.S, idx_basis, rank_info, state_mean)
    elseif basis_mode == "space_local"
        romcfg.basis_centered &&
            error("basis_centered=true is currently supported only for basis_mode=\"global\". Set basis_centered=false for space-local POD.")
        # State ordering is [u; temp], so concatenate coordinates accordingly.
        requested_total_rank = min(requested_rank, Kfull)
        if requested_total_rank < requested_rank
            @warn "Requested npod = $(requested_rank), but only $Kfull snapshots are available. Using rfit = $requested_total_rank."
        end
        x_u, y_u = velocity_xycoords(setup; T = Float64)
        x_t, y_t = scalar_xycoords(setup; T = Float64)
        xcoords = vcat(x_u, x_t)
        ycoords = vcat(y_u, y_t)

        part_method = lowercase(romcfg.basis_partition_method)
        nx = romcfg.basis_partition_nx
        ny = romcfg.basis_partition_ny
        nx >= 1 || error("basis_partition_nx must be >= 1, got $nx")
        ny >= 1 || error("basis_partition_ny must be >= 1, got $ny")

        blocks, partinfo = if part_method == "uniform_x"
            xcuts = uniform_cuts(xcoords, nx)
            blocks_from_xycoords(xcoords, nothing; xcuts = xcuts)
        elseif part_method == "uniform_xy"
            xcuts = uniform_cuts(xcoords, nx)
            ycuts = uniform_cuts(ycoords, ny)
            blocks_from_xycoords(xcoords, ycoords; xcuts = xcuts, ycuts = ycuts)
        elseif part_method == "manual_x"
            xcuts = isnothing(romcfg.basis_partition_xcuts) ? Float64[] : romcfg.basis_partition_xcuts
            blocks_from_xycoords(xcoords, nothing; xcuts = xcuts)
        elseif part_method == "manual_xy"
            xcuts = isnothing(romcfg.basis_partition_xcuts) ? Float64[] : romcfg.basis_partition_xcuts
            ycuts = isnothing(romcfg.basis_partition_ycuts) ? Float64[] : romcfg.basis_partition_ycuts
            blocks_from_xycoords(xcoords, ycoords; xcuts = xcuts, ycuts = ycuts)
        else
            error(
                "Unsupported basis_partition_method: $(romcfg.basis_partition_method). " *
                "Use \"uniform_x\", \"uniform_xy\", \"manual_x\", or \"manual_xy\".",
            )
	        end

        strategy = Symbol(romcfg.basis_local_strategy)
        slpod = build_space_local_basis(
            Xbasis,
            blocks;
            requested_total_rank = requested_total_rank,
            energy_fraction = romcfg.basis_local_energy_fraction,
            local_ranks = romcfg.basis_local_ranks,
            strategy = strategy,
            sparse_basis = romcfg.basis_space_local_sparse,
            row_weights = w,
        )
	        if strategy == :S3 && romcfg.basis_local_ranks === nothing && slpod.rtot != requested_rank
	            @info "Space-local POD (S3) selected an effective rank from the energy threshold" (
	                requested_npod = requested_rank,
	                rfit_effective = slpod.rtot,
	                basis_local_energy_fraction = romcfg.basis_local_energy_fraction,
	            )
	        end

	        rom_basis = slpod.rom_basis
	        projector = slpod.projector
	        singular_values = nothing
	        space_local_info = (;
	            basis_mode = "space_local",
	            basis_local_strategy = String(romcfg.basis_local_strategy),
	            basis_partition_method = String(romcfg.basis_partition_method),
	            basis_partition = partinfo,
	            basis_local_ranks = slpod.ranks,
	            basis_local_energy_fraction = romcfg.basis_local_energy_fraction,
	            basis_space_local_sparse = romcfg.basis_space_local_sparse,
	            rtot = slpod.rtot,
	        )
	        return (; rom_basis, projector, singular_values, idx_basis, space_local_info)
	    else
	        error("Unsupported basis_mode: $(romcfg.basis_mode). Use \"global\" or \"space_local\".")
    end
end

function generate_exact_opinf_snapshots(
    rom_basis::AbstractMatrix,
    basis_projector::AbstractMatrix,
    basis_state_mean,
    setup,
    psolver,
    fomcfg::FOMConfig,
    romcfg::ROMConfig;
    dt_train::Real,
)
    rfit = size(rom_basis, 2)
    T = eltype(rom_basis)
    Atrain = romcfg.include_constant_term ?
        rank_ensuring_states_affine_quadratic(rfit; T = T) :
        rank_ensuring_states_quadratic(rfit; T = T)
    K = size(Atrain, 2)
    Dtrain = zeros(T, rfit, K)

    method = romcfg.exact_opinf_force_explicit_euler ? :explicit_euler : :lmwray3
    method_info = nothing
    state = initial_state(fomcfg, setup, psolver)

    for s in 1:K
        a = @view Atrain[:, s]
        x0 = full_state_from_reduced(rom_basis, a, basis_state_mean)

        # Ensure the training state is physically valid for Dirichlet boundary
        # conditions before taking a one-step difference quotient. This matches
        # the intrusive RHS projection path.
        copy_vector_to_state!(state, x0)
        NS.apply_bc_u!(state.u, 0.0, setup)
        NS.apply_bc_temp!(state.temp, 0.0, setup)
        copy_state_to_vector!(x0, state)

        x1, method_info = one_step_fom_from_state(x0, setup, psolver, fomcfg; dt = dt_train, method = method)
        Dtrain[:, s] = basis_projector * ((x1 - x0) / dt_train)
    end

    return Atrain, Dtrain, method_info
end

function generate_intrusive_rom_snapshots(
    rom_basis::AbstractMatrix,
    basis_projector::AbstractMatrix,
    basis_state_mean,
    setup,
    psolver,
    fomcfg::FOMConfig,
    romcfg::ROMConfig,
)
    rfit = size(rom_basis, 2)
    T = eltype(rom_basis)
    Atrain = romcfg.include_constant_term ?
        rank_ensuring_states_affine_quadratic(rfit; T = T) :
        rank_ensuring_states_quadratic(rfit; T = T)
    K = size(Atrain, 2)
    Dtrain = zeros(T, rfit, K)

    state = initial_state(fomcfg, setup, psolver)
    force = (; u = similar(state.u), temp = similar(state.temp))

    force_cache = NS.get_cache(NS.boussinesq!, setup)
    rk = NS.RKMethods.FE11(; T = T)
    ode_cache = NS.get_cache(rk, state, setup)
    rhs = zeros(T, length(state.u) + length(state.temp))
    params = boussinesq_params(fomcfg)

    for s in 1:K
        a = @view Atrain[:, s]
        x0 = full_state_from_reduced(rom_basis, a, basis_state_mean)

        copy_vector_to_state!(state, x0)
        NS.apply_bc_u!(state.u, 0.0, setup)
        NS.apply_bc_temp!(state.temp, 0.0, setup)

        NS.boussinesq!(force, state, 0.0; setup, cache = force_cache, params...)
        NS.project!(force.u, setup; psolver, ode_cache.p)

        copy_state_to_vector!(rhs, force)
        Dtrain[:, s] = basis_projector * rhs
    end

    method_info = (method = "intrusive_rhs_projection", is_intrusive = true)
    return Atrain, Dtrain, method_info
end

function is_rom_outputs_complete(dir::AbstractString)
    needed = [
        joinpath(dir, "config.toml"),
        joinpath(dir, "fit.jld2"),
        joinpath(dir, "trajectory.jld2"),
        joinpath(dir, "metadata.jld2"),
    ]
    return all(isfile, needed)
end

function _load_local_rom_model(run_id::String, name::AbstractString = run_id)
    fit = load_rom_fit(run_id)
    meta = load_rom_metadata(run_id)
    rom_basis = Matrix(fit.rom_basis)
    basis_projector = Matrix(fit.basis_projector)
    basis_state_mean = if hasproperty(fit, :basis_state_mean) && !isnothing(fit.basis_state_mean)
        Vector(fit.basis_state_mean)
    else
        nothing
    end
    mean_projection = isnothing(basis_state_mean) ? nothing : basis_projector * basis_state_mean
    c = hasproperty(fit, :c) ? fit.c : nothing
    return (;
        run_id,
        name = String(name),
        rfit = size(rom_basis, 2),
        rom_basis,
        basis_projector,
        basis_state_mean,
        mean_projection,
        c = isnothing(c) ? nothing : Vector(c),
        L = Matrix(fit.L),
        Q = Array(fit.Q),
        include_constant_term = !isnothing(c),
        source_training_start_time = hasproperty(meta, :training_start_time) ? meta.training_start_time : nothing,
        source_training_end_time = hasproperty(meta, :training_end_time) ? meta.training_end_time : nothing,
    )
end

function _project_state_to_local(model, x::AbstractVector)
    a = model.basis_projector * x
    isnothing(model.mean_projection) || (a .-= model.mean_projection)
    return a
end

function _reconstruct_local_state!(x::AbstractVector, model, a::AbstractVector)
    length(a) == model.rfit || error("Reduced state length mismatch for local model $(model.name)")
    mul!(x, model.rom_basis, a)
    isnothing(model.basis_state_mean) || (x .+= model.basis_state_mean)
    return x
end

function _weighted_relative_state_difference(xnew::AbstractVector, xold::AbstractVector, w::AbstractVector)
    sdiff = 0.0
    sold = 0.0
    @inbounds for i in eachindex(xnew, xold, w)
        wi = Float64(w[i])
        di = Float64(xnew[i]) - Float64(xold[i])
        oi = Float64(xold[i])
        sdiff += wi * di * di
        sold += wi * oi * oi
    end
    return sqrt(sdiff) / max(sqrt(sold), eps(Float64))
end

function _normalized_lz_from_state_vector!(x::AbstractVector, state, setup, angular_cache, t::Real)
    copy_vector_to_state!(state, x)
    NS.apply_bc_u!(state.u, t, setup)
    return normalized_angular_momentum_z(state.u, setup, angular_cache)
end

function _same_fom_value(a, b)
    if a isa AbstractFloat || b isa AbstractFloat
        return isapprox(Float64(a), Float64(b); rtol = 1e-12, atol = 0.0)
    end
    return a == b
end

function _fom_cross_prediction_mismatches(training_cfg::FOMConfig, prediction_cfg::FOMConfig)
    ignored = Set([:tsim, :nsnapshot])
    mismatches = Symbol[]
    for field in fieldnames(FOMConfig)
        field in ignored && continue
        _same_fom_value(getfield(training_cfg, field), getfield(prediction_cfg, field)) ||
            push!(mismatches, field)
    end
    return mismatches
end

function _check_local_rom_parent_compatible!(
    model,
    local_meta,
    prediction_parent_fom_id::String,
    prediction_fomcfg::FOMConfig,
    switchcfg::SwitchedROMConfig,
)
    hasproperty(local_meta, :parent_fom_id) || return nothing
    training_parent_fom_id = String(local_meta.parent_fom_id)
    training_parent_fom_id == prediction_parent_fom_id && return nothing

    if !switchcfg.allow_parent_fom_mismatch
        error("Local ROM $(model.run_id) was trained from $(training_parent_fom_id), not $(prediction_parent_fom_id)")
    end

    training_fomcfg = load_fom_config_from_run(training_parent_fom_id)
    mismatches = _fom_cross_prediction_mismatches(training_fomcfg, prediction_fomcfg)
    if !isempty(mismatches)
        error(
            "Local ROM $(model.run_id) was trained from $(training_parent_fom_id), " *
            "but the prediction FOM $(prediction_parent_fom_id) differs in $(join(string.(mismatches), ", ")). " *
            "Cross-parent switched ROM prediction is only allowed when FOM configs match except for tsim/nsnapshot.",
        )
    end

    @info "Allowing switched-ROM local model trained on compatible shorter/longer FOM" (
        model = model.name,
        local_rom_id = model.run_id,
        training_parent_fom_id,
        prediction_parent_fom_id,
    )
    return nothing
end

function _switched_target_index(active::Integer, lz_normalized::Real, t::Real, cfg::SwitchedROMConfig)
    if active == 1
        return t >= cfg.pre_to_transition_time ? 2 : 1
    elseif active == 2
        return t >= cfg.transition_to_post_time ? 3 : 2
    elseif active == 3
        return 3
    end
    error("Unsupported active switched ROM index: $active")
end

function _advance_local_rom_interval!(
    a::AbstractVector,
    model,
    Δt::Real;
    safety::Real,
    growthfac::Real,
    dtmin::Real = 1e-10,
    absolute_time::Real = 0.0,
)
    Δt_float = Float64(Δt)
    Δt_float >= 0.0 || error("Cannot advance a local ROM by a negative time interval")
    Δt_float == 0.0 && return a

    r = length(a)
    T = eltype(a)
    atrial = similar(a)
    k1 = similar(a)
    k2 = similar(a)
    k3 = similar(a)
    k4 = similar(a)
    atmp = similar(a)
    Jtmp = Matrix{T}(undef, r, r)

    tlocal = 0.0
    isubstep = 0
    while tlocal < prevfloat(Δt_float)
        dtcap = Δt_float - tlocal
        dt = min(propose_rom_timestep(a, model.L, model.Q; dtmax = dtcap, safety = safety), dtcap)

        accepted = false
        while !accepted
            copyto!(atrial, a)
            if model.include_constant_term
                rk4_step!(atrial, dt, model.c, model.L, model.Q, k1, k2, k3, k4, atmp)
            else
                rk4_step!(atrial, dt, model.L, model.Q, k1, k2, k3, k4, atmp)
            end

            bad = !all(isfinite, atrial)
            na = norm(a)
            ntrial = norm(atrial)
            growth_bad = ntrial > growthfac * max(na, 1e-12)

            if bad || growth_bad
                rom_jacobian!(Jtmp, a, model.L, model.Q)
                @info "Rejected switched-ROM local step" (
                    model = model.name,
                    t = Float64(absolute_time) + tlocal,
                    dt = dt,
                    jacnorm = opnorm(Jtmp, 1),
                    norma = na,
                    normtrial = ntrial,
                    bad = bad,
                    growth_bad = growth_bad,
                )

                dt *= 0.5
                if dt < dtmin
                    error("Switched ROM timestep became too small at t = $(Float64(absolute_time) + tlocal)")
                end
            else
                copyto!(a, atrial)
                tlocal += dt
                isubstep += 1
                accepted = true
            end
        end
    end

    return a
end

function is_switched_rom_outputs_complete(dir::AbstractString)
    return is_rom_outputs_complete(dir) && isfile(joinpath(dir, "switched_rom_info.txt"))
end

function run_switched_rom(parent_fom_id::String, switchcfg::SwitchedROMConfig; force::Bool = false, skip_if_exists::Bool = true)
    validate_switched_rom_config(switchcfg)
    run_id = switched_rom_id(parent_fom_id, switchcfg)
    dir = rom_dir(run_id)
    run_start_ns = time_ns()

    if isdir(dir)
        if skip_if_exists
            if is_switched_rom_outputs_complete(dir)
                @info "Switched ROM run already exists, reusing $run_id"
                return run_id
            else
                @warn "Switched ROM run directory exists but outputs missing; rerunning" run_id dir
                rm(dir; recursive = true, force = true)
            end
        elseif !force
            error("Switched ROM directory already exists: $dir")
        else
            rm(dir; recursive = true, force = true)
        end
    end

    ensure_run_dirs!(dir)
    save_switched_rom_config(joinpath(dir, "config.toml"), switchcfg)
    open(joinpath(dir, "parent_fom_id.txt"), "w") do io
        write(io, parent_fom_id)
        write(io, "\n")
    end

    fomcfg = load_fom_config_from_run(parent_fom_id)
    snapshots = load_fom_snapshots(parent_fom_id)
    Kfull = size(snapshots.state, 2)
    prediction_window = rom_training_window(fomcfg, Kfull, switchcfg.start_time, switchcfg.end_time)
    prediction_indices = prediction_window.snapshot_indices
    prediction_times = collect(prediction_window.times)
    K = length(prediction_indices)
    K >= 2 || error("Switched ROM prediction window must include at least two snapshots")

    models = [
        _load_local_rom_model(run_id, switchcfg.regime_names[i])
        for (i, run_id) in enumerate(switchcfg.local_rom_ids)
    ]
    for model in models
        local_meta = load_rom_metadata(model.run_id)
        _check_local_rom_parent_compatible!(model, local_meta, parent_fom_id, fomcfg, switchcfg)
    end

    initial_active = findfirst(==(switchcfg.initial_regime), switchcfg.regime_names)
    isnothing(initial_active) && error("Could not find initial_regime=$(switchcfg.initial_regime)")
    active = Int(initial_active)

    setup = getsetup(fomcfg; backend = KernelAbstractions.CPU())
    psolver = NS.default_psolver(setup)
    state = initial_state(fomcfg, setup, psolver)
    angular_cache = angular_momentum_z_cache(setup, state.u, fomcfg)
    state_weights = state_mass_vector(setup, state.u, state.temp)

    rmax = maximum(model.rfit for model in models)
    a_rom = fill(NaN, rmax, K)
    a_exact = fill(NaN, rmax, K)
    active_model_index = zeros(Int, K)
    switch_metric = zeros(Float64, K)

    xbuf = Vector{Float64}(undef, size(snapshots.state, 1))
    xproj = similar(xbuf)
    first_state = @view snapshots.state[:, first(prediction_indices)]
    a = _project_state_to_local(models[active], first_state)

    switch_events = NamedTuple[]

    function record_snapshot!(k::Integer, active_idx::Integer, ak::AbstractVector)
        model = models[active_idx]
        r = model.rfit
        active_model_index[k] = active_idx
        copyto!(@view(a_rom[1:r, k]), ak)
        a_fom = _project_state_to_local(model, @view(snapshots.state[:, prediction_indices[k]]))
        copyto!(@view(a_exact[1:r, k]), a_fom)
        _reconstruct_local_state!(xbuf, model, ak)
        switch_metric[k] = _normalized_lz_from_state_vector!(xbuf, state, setup, angular_cache, prediction_times[k])
        return nothing
    end

    advance_runtime_seconds = 0.0
    online_runtime_seconds = @elapsed begin
        record_snapshot!(1, active, a)

        prog = Progress(K - 1; desc = "Running switched ROM")
        for k in 1:(K - 1)
            t0 = prediction_times[k]
            t1 = prediction_times[k + 1]
            model = models[active]
            advance_runtime_seconds += @elapsed _advance_local_rom_interval!(
                a,
                model,
                t1 - t0;
                safety = switchcfg.safety,
                growthfac = switchcfg.growthfac,
                absolute_time = t0,
            )

            _reconstruct_local_state!(xbuf, model, a)
            lz_before = _normalized_lz_from_state_vector!(xbuf, state, setup, angular_cache, t1)
            target = _switched_target_index(active, lz_before, t1, switchcfg)
            if target != active
                old_active = active
                new_model = models[target]
                anew = _project_state_to_local(new_model, xbuf)
                _reconstruct_local_state!(xproj, new_model, anew)
                lz_after = _normalized_lz_from_state_vector!(xproj, state, setup, angular_cache, t1)
                jump_rel = _weighted_relative_state_difference(xproj, xbuf, state_weights)
                push!(
                    switch_events,
                    (
                        time = t1,
                        from_index = old_active,
                        to_index = target,
                        from_regime = switchcfg.regime_names[old_active],
                        to_regime = switchcfg.regime_names[target],
                        lz_before_projection = lz_before,
                        lz_after_projection = lz_after,
                        projection_jump_rel = jump_rel,
                    ),
                )
                @info "Switched ROM regime" switch_events[end]
                active = target
                a = anew
            end

            record_snapshot!(k + 1, active, a)
            next!(
                prog;
                showvalues = [
                    ("time", t1),
                    ("regime", switchcfg.regime_names[active]),
                    ("Lz_norm", switch_metric[k + 1]),
                    ("switches", length(switch_events)),
                ],
            )
        end
    end

    fit_data = (
        switched = true,
        switching_variable = "angular_momentum_z_normalized",
        local_rom_ids = switchcfg.local_rom_ids,
        local_model_names = switchcfg.regime_names,
        local_rfits = [model.rfit for model in models],
    )
    traj_data = (
        a_exact = a_exact,
        dadt_exact = nothing,
        a_rom = a_rom,
        times = prediction_times,
        fom_snapshot_indices = collect(prediction_indices),
        active_model_index = active_model_index,
        switch_metric = switch_metric,
        switch_events = switch_events,
    )
    thresholds = (
        switching_policy = switchcfg.switching_policy,
        min_switch_time = switchcfg.min_switch_time,
        pre_to_transition_threshold = switchcfg.pre_to_transition_threshold,
        transition_to_post_threshold = switchcfg.transition_to_post_threshold,
        transition_to_pre_threshold = switchcfg.transition_to_pre_threshold,
        post_to_transition_threshold = switchcfg.post_to_transition_threshold,
        pre_to_transition_time = switchcfg.pre_to_transition_time,
        transition_to_post_time = switchcfg.transition_to_post_time,
    )
    local_training_end_times = [
        model.source_training_end_time for model in models
        if !isnothing(model.source_training_end_time)
    ]
    local_training_end_time = isempty(local_training_end_times) ? nothing : maximum(Float64.(local_training_end_times))
    run_setup_and_online_seconds = (time_ns() - run_start_ns) / 1e9
    metadata = (
        created_at = string(Dates.now()),
        runtime_seconds = online_runtime_seconds,
        online_runtime_seconds = online_runtime_seconds,
        advance_runtime_seconds = advance_runtime_seconds,
        setup_and_online_runtime_seconds = run_setup_and_online_seconds,
        timing_scope = "online_runtime_seconds includes switched ROM prediction, switching, reconstruction, and FOM projection diagnostics; advance_runtime_seconds includes only reduced-coordinate RK advances",
        parent_fom_id = parent_fom_id,
        run_tag = switchcfg.run_tag,
        switched_rom = true,
        switching_variable = "angular_momentum_z_normalized",
        switching_policy = switchcfg.switching_policy,
        switch_thresholds = thresholds,
        allow_backward_switches = switchcfg.allow_backward_switches,
        pre_to_transition_time = switchcfg.pre_to_transition_time,
        transition_to_post_time = switchcfg.transition_to_post_time,
        allow_parent_fom_mismatch = switchcfg.allow_parent_fom_mismatch,
        local_rom_ids = switchcfg.local_rom_ids,
        local_model_names = switchcfg.regime_names,
        local_rfits = [model.rfit for model in models],
        local_training_end_time = local_training_end_time,
        rfit = rmax,
        initial_regime = switchcfg.initial_regime,
        training_start_time_requested = prediction_window.requested_start_time,
        training_end_time_requested = prediction_window.requested_end_time,
        training_start_time = prediction_window.actual_start_time,
        training_end_time = prediction_window.end_time,
        training_duration = prediction_window.duration,
        nswitches = length(switch_events),
        switch_events = switch_events,
    )

    open(joinpath(dir, "rom_dim.txt"), "w") do io
        write(io, string(rmax))
        write(io, "\n")
    end
    open(joinpath(dir, "switched_rom_info.txt"), "w") do io
        write(io, "switching_variable = angular_momentum_z_normalized\n")
        write(io, "start_time_requested = $(prediction_window.requested_start_time)\n")
        write(io, "end_time_requested = $(prediction_window.requested_end_time)\n")
        write(io, "start_time = $(prediction_window.actual_start_time)\n")
        write(io, "end_time = $(prediction_window.end_time)\n")
        write(io, "switching_policy = $(switchcfg.switching_policy)\n")
        write(io, "allow_parent_fom_mismatch = $(switchcfg.allow_parent_fom_mismatch)\n")
        write(io, "allow_backward_switches = $(switchcfg.allow_backward_switches)\n")
        write(io, "pre_to_transition_time = $(switchcfg.pre_to_transition_time)\n")
        write(io, "transition_to_post_time = $(switchcfg.transition_to_post_time)\n")
        write(io, "local_training_end_time = $(local_training_end_time)\n")
        write(io, "runtime_seconds = $(online_runtime_seconds)\n")
        write(io, "advance_runtime_seconds = $(advance_runtime_seconds)\n")
        write(io, "setup_and_online_runtime_seconds = $(run_setup_and_online_seconds)\n")
        write(io, "rmax = $rmax\n")
        write(io, "nswitches = $(length(switch_events))\n")
        for (i, model) in enumerate(models)
            write(io, "model_$i = $(model.name), $(model.run_id), r=$(model.rfit), train=$(model.source_training_start_time)..$(model.source_training_end_time)\n")
        end
        write(io, "threshold_pre_to_transition = $(switchcfg.pre_to_transition_threshold)\n")
        write(io, "threshold_transition_to_post = $(switchcfg.transition_to_post_threshold)\n")
        for ev in switch_events
            write(io, "switch = t=$(ev.time), $(ev.from_regime)->$(ev.to_regime), Lz=$(ev.lz_before_projection), jump=$(ev.projection_jump_rel)\n")
        end
    end

    save_object(joinpath(dir, "fit.jld2"), fit_data)
    save_object(joinpath(dir, "trajectory.jld2"), traj_data)
    save_object(joinpath(dir, "metadata.jld2"), metadata)

    if isdefined(@__MODULE__, :register_switched_rom!)
        register_switched_rom!(run_id, parent_fom_id, switchcfg, fomcfg)
    end
    @info "Saved switched ROM run to $dir"

    return run_id
end

function run_rom(parent_fom_id::String, romcfg::ROMConfig; force::Bool = false, skip_if_exists::Bool = true)
    run_id = rom_id(parent_fom_id, romcfg)
    dir = rom_dir(run_id)

    if isdir(dir)
        if skip_if_exists
            if is_rom_outputs_complete(dir)
                open(joinpath(dir, "parent_fom_id.txt"), "w") do io
                    write(io, parent_fom_id)
                    write(io, "\n")
                end
                if !isfile(joinpath(dir, "rom_dim.txt")) && isfile(joinpath(dir, "metadata.jld2"))
                    meta = load_object(joinpath(dir, "metadata.jld2"))
                    if hasproperty(meta, :rfit)
                        open(joinpath(dir, "rom_dim.txt"), "w") do io
                            write(io, string(meta.rfit))
                            write(io, "\n")
                        end
                    end
                end
                if !isfile(joinpath(dir, "rom_rank_info.txt")) && isfile(joinpath(dir, "metadata.jld2"))
                    meta = load_object(joinpath(dir, "metadata.jld2"))
                    # Best-effort reconstruction for older runs.
                    basis_mode = hasproperty(meta, :basis_mode) ? string(meta.basis_mode) : "unknown"
                    rank_selection = hasproperty(meta, :rank_selection) ? string(meta.rank_selection) : "unknown"
                    requested_npod = hasproperty(meta, :requested_npod) ? meta.requested_npod : nothing
                    energy_fraction = hasproperty(meta, :energy_fraction) ? meta.energy_fraction : nothing
                    energy_retained =
                        hasproperty(meta, :rank_info) && meta.rank_info !== nothing && hasproperty(meta.rank_info, :energy_retained) ?
                        meta.rank_info.energy_retained : nothing
                    rfit = hasproperty(meta, :rfit) ? meta.rfit : nothing
                    open(joinpath(dir, "rom_rank_info.txt"), "w") do io
                        write(io, "basis_mode = $basis_mode\n")
                        write(io, "rank_selection = $rank_selection\n")
                        if requested_npod !== nothing
                            write(io, "requested_npod = $requested_npod\n")
                        end
                        if energy_fraction !== nothing
                            write(io, "energy_fraction = $energy_fraction\n")
                        end
                        if energy_retained !== nothing
                            write(io, "energy_retained = $(round(energy_retained; digits = 10))\n")
                        end
                        if rfit !== nothing
                            write(io, "rfit = $rfit\n")
                        end
                    end
                end
                @info "ROM run already exists, reusing $run_id"
                return run_id
            else
                @warn "ROM run directory exists but outputs missing; rerunning" run_id dir
                rm(dir; recursive = true, force = true)
            end
        elseif !force
            error("ROM directory already exists: $dir")
        else
            rm(dir; recursive = true, force = true)
        end
    end

    ensure_run_dirs!(dir)
    save_rom_config(joinpath(dir, "config.toml"), romcfg)
    open(joinpath(dir, "parent_fom_id.txt"), "w") do io
        write(io, parent_fom_id)
        write(io, "\n")
    end

    fomcfg = load_fom_config_from_run(parent_fom_id)
    snapshots = load_fom_snapshots(parent_fom_id)

    Kfull = size(snapshots.state, 2)
    train_window = rom_training_window(fomcfg, Kfull, romcfg.training_start_time, romcfg.training_end_time)
    train_indices = train_window.snapshot_indices
    snapshots_state_train = @view snapshots.state[:, train_indices]
    snapshots_dstate_train = @view snapshots.dstate[:, train_indices]
    @info "Selected ROM training window" (
        requested_start_time = train_window.requested_start_time,
        requested_end_time = train_window.requested_end_time,
        actual_start_time = train_window.actual_start_time,
        end_time = train_window.end_time,
        nsnapshots = length(train_indices),
        first_snapshot = first(train_indices),
        last_snapshot = last(train_indices),
    )

    use_intrusive_rom_effective = romcfg.use_intrusive_rom
    use_exact_opinf_effective = !use_intrusive_rom_effective && romcfg.use_exact_opinf
    if use_exact_opinf_effective && !romcfg.exact_opinf_force_explicit_euler && romcfg.fallback_to_standard_opinf
        @warn "Falling back to trajectory-based OpInf because explicit Euler was not forced." (
            exact_opinf_force_explicit_euler = romcfg.exact_opinf_force_explicit_euler,
            fallback_to_standard_opinf = romcfg.fallback_to_standard_opinf,
        )
        use_exact_opinf_effective = false
    end
    train_mode = use_intrusive_rom_effective ? :intrusive : (use_exact_opinf_effective ? :exact_opinf : :traj_opinf)

    Atrain = nothing
    Dtrain = nothing
    method_info = nothing
    train_dt = nothing

    # Build ROM basis (always on CPU for deterministic coordinate mapping).
    setup_basis = getsetup(fomcfg; backend = KernelAbstractions.CPU())
    basis = build_rom_basis(
        snapshots_state_train,
        setup_basis,
        romcfg;
        requested_rank = romcfg.npod,
        snapshot_indices = train_indices,
    )
    rom_basis = basis.rom_basis
    basis_projector = basis.projector
    basis_state_mean = hasproperty(basis, :state_mean) ? basis.state_mean : nothing
    rfit = size(rom_basis, 2)
    basis_space_local_info = hasproperty(basis, :space_local_info) ? basis.space_local_info : nothing
    rank_info = hasproperty(basis, :rank_info) ? basis.rank_info : nothing
    a_exact = project_snapshots_with_mean(basis_projector, snapshots_state_train, basis_state_mean)
    dadt_exact = basis_projector * snapshots_dstate_train

    if !isnothing(rank_info) && hasproperty(rank_info, :energy_retained)
        @info "ROM POD energy retained" (
            rank_selection = rank_info.rank_selection,
            energy_basis = hasproperty(rank_info, :energy_basis) ? rank_info.energy_basis : "unknown",
            requested_npod = rank_info.requested_npod,
            rfit = rank_info.rfit,
            energy_retained = rank_info.energy_retained,
        )
    end

    if train_mode in (:intrusive, :exact_opinf)
        backend = getbackend()
        ensure_cuda_psolver!(backend)
        setup = getsetup(fomcfg; backend = backend)
        psolver = rayleigh_psolver(fomcfg, setup)

        if train_mode == :intrusive
            Atrain, Dtrain, method_info = generate_intrusive_rom_snapshots(
                rom_basis,
                basis_projector,
                basis_state_mean,
                setup,
                psolver,
                fomcfg,
                romcfg,
            )
            fit = romcfg.include_constant_term ?
                fit_affine_quadratic_rom_exact_symmetric(Atrain, Dtrain) :
                fit_quadratic_rom_exact_symmetric(Atrain, Dtrain)
            @info "Intrusive ROM training diagnostics" (
                K_train = size(Atrain, 2),
                nf = romcfg.include_constant_term ? n_features_affine_quadratic(rfit) : (rfit + n_unique_quadratic(rfit)),
                condPhi = fit.σ[1] / max(fit.σ[end], eps(eltype(fit.σ))),
                resnorm = fit.resnorm,
                relres = fit.relres,
                method = method_info,
            )
        else
            snapshot_dt = fomcfg.tsim / fomcfg.nsnapshot
            estimated_dt = estimate_exact_opinf_dt_matlab(a_exact, dadt_exact; fallback_dt = snapshot_dt)
            train_dt = min(estimated_dt, snapshot_dt)
            if train_dt < estimated_dt
                @info "Capped exact-OpInf training step at the FOM snapshot spacing" estimated_dt train_dt snapshot_dt
            else
                @info "Estimated exact-OpInf training step" train_dt snapshot_dt
            end
            Atrain, Dtrain, method_info = generate_exact_opinf_snapshots(
                rom_basis,
                basis_projector,
                basis_state_mean,
                setup,
                psolver,
                fomcfg,
                romcfg;
                dt_train = train_dt,
            )

            fit = romcfg.include_constant_term ?
                fit_affine_quadratic_rom_exact_symmetric(Atrain, Dtrain) :
                fit_quadratic_rom_exact_symmetric(Atrain, Dtrain)
            if !romcfg.exact_opinf_force_explicit_euler
                @warn "Exact OpInf training did not use explicit Euler; exact recovery assumptions are not fully satisfied." method = method_info
            end

            @info "Exact OpInf training diagnostics" (
                K_train = size(Atrain, 2),
                nf = romcfg.include_constant_term ? n_features_affine_quadratic(rfit) : (rfit + n_unique_quadratic(rfit)),
                condPhi = fit.σ[1] / max(fit.σ[end], eps(eltype(fit.σ))),
                resnorm = fit.resnorm,
                relres = fit.relres,
                method = method_info,
            )
        end
    else
        K = length(train_indices)
        rmax = romcfg.include_constant_term ?
            max_affine_quadratic_rank_symmetric(K) :
            max_quadratic_rank_symmetric(K)
        rdesired = size(a_exact, 1)
        rfit_traj = min(rdesired, rmax)
        if rfit_traj < rdesired
            @warn "Selected basis rank = $rdesired, but with K = $K only up to r = $rmax is identifiable for trajectory OpInf. Using rfit = $rfit_traj."
        end

        a_train = @view a_exact[1:rfit_traj, :]
        dadt_train = @view dadt_exact[1:rfit_traj, :]
        fit = romcfg.include_constant_term ?
            fit_affine_quadratic_rom_symmetric(a_train, dadt_train; lambda = romcfg.lambda) :
            fit_quadratic_rom_symmetric(a_train, dadt_train; lambda = romcfg.lambda)

        # Truncate basis and coordinates consistently for rollout.
        rom_basis = @view rom_basis[:, 1:rfit_traj]
        basis_projector = @view basis_projector[1:rfit_traj, :]
        a_exact = @view a_exact[1:rfit_traj, :]
        dadt_exact = @view dadt_exact[1:rfit_traj, :]
        rfit = rfit_traj
    end

    a0 = copy(a_exact[:, 1])
    rollout_duration = train_window.duration
    rollout_nsave = length(train_indices) - 1
    a_rom = if romcfg.include_constant_term
        simulate_rom_adaptive(
            a0,
            fit.c,
            fit.L,
            fit.Q;
            tsim = rollout_duration,
            nsave = rollout_nsave,
            safety = romcfg.safety,
            growthfac = romcfg.growthfac,
        )
    else
        simulate_rom_adaptive(
            a0,
            fit.L,
            fit.Q;
            tsim = rollout_duration,
            nsave = rollout_nsave,
            safety = romcfg.safety,
            growthfac = romcfg.growthfac,
        )
    end

    condphi = fit.σ[1] / max(fit.σ[end], eps(eltype(fit.σ)))
    diagnostics = (
        rfit = rfit,
        condphi = condphi,
        maxabsL = maximum(abs, fit.L),
        maxabsQ = maximum(abs, fit.Q),
    )

    fit_data = (
        rom_basis = rom_basis,
        basis_projector = basis_projector,
        basis_state_mean = basis_state_mean,
        basis_idx = basis.idx_basis,
        c = romcfg.include_constant_term ? fit.c : nothing,
        L = fit.L,
        Q = fit.Q,
        B = fit.B,
        σ = fit.σ,
        pairs = fit.pairs,
        train_mode = String(train_mode),
        use_intrusive_rom = use_intrusive_rom_effective,
        use_intrusive_rom_config = romcfg.use_intrusive_rom,
        use_exact_opinf = use_exact_opinf_effective,
        use_exact_opinf_config = romcfg.use_exact_opinf,
        train_dt = use_exact_opinf_effective ? train_dt : nothing,
        train_residual_norm = hasproperty(fit, :resnorm) ? fit.resnorm : nothing,
        train_residual_rel = hasproperty(fit, :relres) ? fit.relres : nothing,
        Atrain = Atrain,
        Dtrain = Dtrain,
        train_method = method_info,
        training_start_time_requested = train_window.requested_start_time,
        training_end_time_requested = train_window.requested_end_time,
        training_start_time = train_window.actual_start_time,
        training_end_time = train_window.end_time,
        diagnostics = diagnostics,
    )
    traj_data = (
        a_exact = a_exact,
        dadt_exact = dadt_exact,
        a_rom = a_rom,
        times = train_window.times,
        fom_snapshot_indices = collect(train_indices),
    )
    metadata = (
        created_at = string(Dates.now()),
        parent_fom_id = parent_fom_id,
        run_tag = romcfg.run_tag,
        requested_npod = romcfg.npod,
        rfit = rfit,
        basis_mode = lowercase(romcfg.basis_mode),
        rank_selection = isnothing(rank_info) ? "unknown" : rank_info.rank_selection,
        rank_info = rank_info,
        basis_centered = !isnothing(basis_state_mean),
        basis_space_local_info = basis_space_local_info,
        diagnostics = diagnostics,
        train_mode = String(train_mode),
        use_intrusive_rom = use_intrusive_rom_effective,
        use_intrusive_rom_config = romcfg.use_intrusive_rom,
        use_exact_opinf = use_exact_opinf_effective,
        use_exact_opinf_config = romcfg.use_exact_opinf,
        exact_opinf_force_explicit_euler = romcfg.exact_opinf_force_explicit_euler,
        exact_opinf_assumptions_satisfied = use_exact_opinf_effective && romcfg.exact_opinf_force_explicit_euler,
        include_constant_term = romcfg.include_constant_term,
        training_start_time_requested = train_window.requested_start_time,
        training_end_time_requested = train_window.requested_end_time,
        training_start_time = train_window.actual_start_time,
        training_end_time = train_window.end_time,
        training_duration = train_window.duration,
        dt_train = use_exact_opinf_effective ? train_dt : nothing,
        K_train = train_mode in (:intrusive, :exact_opinf) ?
            (romcfg.include_constant_term ? n_features_affine_quadratic(rfit) : (rfit + n_unique_quadratic(rfit))) :
            size(a_exact, 2),
        basis_snapshot_start = first(basis.idx_basis),
        basis_snapshot_end = last(basis.idx_basis),
        basis_nsnapshot = length(basis.idx_basis),
    )

    open(joinpath(dir, "rom_dim.txt"), "w") do io
        # Effective ROM dimension (can differ from requested npod for space-local S3).
        write(io, string(rfit))
        write(io, "\n")
    end
    open(joinpath(dir, "rom_rank_info.txt"), "w") do io
        write(io, "basis_mode = $(lowercase(romcfg.basis_mode))\n")
        if !isnothing(rank_info)
            write(io, "rank_selection = $(rank_info.rank_selection)\n")
            if hasproperty(rank_info, :energy_basis)
                write(io, "energy_basis = $(rank_info.energy_basis)\n")
            end
        else
            write(io, "rank_selection = unknown\n")
        end
        write(io, "requested_npod = $(romcfg.npod)\n")
        write(io, "energy_fraction = $(romcfg.energy_fraction)\n")
        write(io, "training_start_time_requested = $(train_window.requested_start_time)\n")
        write(io, "training_end_time_requested = $(train_window.requested_end_time)\n")
        write(io, "training_start_time = $(train_window.actual_start_time)\n")
        write(io, "training_end_time = $(train_window.end_time)\n")
        if !isnothing(rank_info) && hasproperty(rank_info, :energy_retained)
            write(io, "energy_retained = $(round(rank_info.energy_retained; digits = 10))\n")
        end
        if !isnothing(rank_info) && hasproperty(rank_info, :raw_energy_retained) && !isnothing(rank_info.raw_energy_retained)
            write(io, "raw_energy_retained = $(round(rank_info.raw_energy_retained; digits = 10))\n")
        end
        if !isnothing(rank_info) && hasproperty(rank_info, :fluctuation_energy_retained) && !isnothing(rank_info.fluctuation_energy_retained)
            write(io, "fluctuation_energy_retained = $(round(rank_info.fluctuation_energy_retained; digits = 10))\n")
        end
        write(io, "basis_centered = $(!isnothing(basis_state_mean))\n")
        write(io, "rfit = $(rfit)\n")
    end

    save_object(joinpath(dir, "fit.jld2"), fit_data)
    save_object(joinpath(dir, "trajectory.jld2"), traj_data)
    save_object(joinpath(dir, "metadata.jld2"), metadata)

    register_rom!(run_id, parent_fom_id, romcfg, fomcfg; dt_train = train_dt)
    @info "Saved ROM run to $dir"

    return run_id
end

load_rom_fit(run_id::String) = load_object(joinpath(rom_dir(run_id), "fit.jld2"))
load_rom_trajectory(run_id::String) = load_object(joinpath(rom_dir(run_id), "trajectory.jld2"))
load_rom_metadata(run_id::String) = load_object(joinpath(rom_dir(run_id), "metadata.jld2"))
load_rom_config_from_run(run_id::String) = load_rom_config(joinpath(rom_dir(run_id), "config.toml"))
