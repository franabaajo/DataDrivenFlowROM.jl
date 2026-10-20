using SHA

repo_root() = normpath(joinpath(@__DIR__, "..", ".."))
examples_root() = repo_root() # Backward-compatible name used by the registry code.
dipole_case_name() = "DipoleWallCollision"
default_output_root() = joinpath(repo_root(), "output", dipole_case_name(), "dipolepod")
legacy_output_root() = joinpath(repo_root(), "output", "dipolepod")

function output_root()
    env_root = get(ENV, "DIPOLE_OUTPUT_ROOT", "")
    if !isempty(env_root)
        return mkpath(normpath(env_root))
    end

    preferred = default_output_root()
    if isdir(preferred)
        return mkpath(preferred)
    end

    legacy = legacy_output_root()
    if isdir(legacy)
        @info "Using legacy dipole output root" legacy
        return mkpath(legacy)
    end

    return mkpath(preferred)
end

fom_root() = mkpath(joinpath(output_root(), "fom"))
rom_root() = mkpath(joinpath(output_root(), "rom"))

fom_dir(fom_id::String) = joinpath(fom_root(), fom_id)
rom_dir(rom_id::String) = joinpath(rom_root(), rom_id)

function ensure_run_dirs!(dir::AbstractString)
    mkpath(dir)
    mkpath(joinpath(dir, "plots"))
    mkpath(joinpath(dir, "movies"))
    return dir
end

function fmt_token(x)
    if isnothing(x)
        return "uniform"
    end
    s = lowercase(string(x))
    s = replace(s, "." => "p", "-" => "m", "+" => "")
    return s
end

function tag_token(x)
    s = fmt_token(x)
    s = replace(s, r"[^a-z0-9_]+" => "_")
    s = strip(s, '_')
    return isempty(s) ? "tagged" : s
end

short_hash(s::AbstractString) = bytes2hex(sha1(s))[1:8]

function parent_fom_hash_token(parent_fom_id::AbstractString)
    m = match(r"_([0-9a-f]{8})$", parent_fom_id)
    if !isnothing(m)
        return m.captures[1]
    end
    @warn "Could not parse parent FOM hash token; falling back to hashed parent id" parent_fom_id
    return short_hash(parent_fom_id)
end

function fom_signature(cfg::FOMConfig; tag::Union{Nothing,AbstractString} = nothing)
    parts = [
        "viscosity=$(cfg.viscosity)",
        "r0=$(cfg.r0)",
        "p1=$(cfg.p1)",
        "p2=$(cfg.p2)",
        "initialenergy=$(cfg.initialenergy)",
        "n=$(cfg.n)",
        "stretch=$(cfg.stretch)",
        "lims=$(cfg.lims)",
        "tsim=$(cfg.tsim)",
        "cfl=$(cfg.cfl)",
        "nsnapshot=$(cfg.nsnapshot)",
    ]
    if !isnothing(tag)
        push!(parts, "tag=$(tag)")
    end
    join(parts, "|")
end

function rom_signature(parent_fom_id::String, cfg::ROMConfig)
    parts = String[
        "parent=$(parent_fom_id)",
        "npod=$(cfg.npod)",
        "use_intrusive_rom=$(cfg.use_intrusive_rom)",
        "use_exact_opinf=$(cfg.use_exact_opinf)",
        "force_explicit_euler=$(cfg.exact_opinf_force_explicit_euler)",
        "fallback_to_standard_opinf=$(cfg.fallback_to_standard_opinf)",
        "safety=$(cfg.safety)",
        "growthfac=$(cfg.growthfac)",
    ]
    # Lambda is only used for trajectory-based operator inference.
    method = rom_method_token(cfg)
    if method == "trajopinf"
        push!(parts, "lambda=$(cfg.lambda)")
    end
    rank_sel = lowercase(strip(cfg.rank_selection))
    if rank_sel in ("energy_fraction", "energyfraction", "energy", "fraction", "eta")
        push!(parts, "rank_selection=$(rank_sel)")
        push!(parts, "energy_fraction=$(cfg.energy_fraction)")
    end
    basis_mode = lowercase(cfg.basis_mode)
    push!(parts, "basis_mode=$(basis_mode)")
    if basis_mode == "space_local"
        push!(parts, "basis_local_strategy=$(lowercase(cfg.basis_local_strategy))")
        push!(parts, "basis_local_energy_fraction=$(cfg.basis_local_energy_fraction)")
        push!(parts, "basis_local_ranks=$(cfg.basis_local_ranks)")
        push!(parts, "basis_partition_method=$(lowercase(cfg.basis_partition_method))")
        push!(parts, "basis_partition_nx=$(cfg.basis_partition_nx)")
        push!(parts, "basis_partition_ny=$(cfg.basis_partition_ny)")
        push!(parts, "basis_partition_xcuts=$(cfg.basis_partition_xcuts)")
        push!(parts, "basis_partition_ycuts=$(cfg.basis_partition_ycuts)")
        push!(parts, "basis_space_local_sparse=$(cfg.basis_space_local_sparse)")
    end
    if !isnothing(cfg.run_tag)
        push!(parts, "run_tag=$(cfg.run_tag)")
    end
    return join(parts, "|")
end

function rom_method_token(cfg::ROMConfig)
    if cfg.use_intrusive_rom
        return "intrusive"
    end
    if cfg.use_exact_opinf
        # If explicit Euler is not forced and fallback is enabled, the workflow switches
        # to trajectory-based OpInf.
        if !cfg.exact_opinf_force_explicit_euler && cfg.fallback_to_standard_opinf
            return "trajopinf"
        end
        return "exactopinf"
    end
    return "trajopinf"
end

function fom_id(
    cfg::FOMConfig;
    tag::Union{Nothing,AbstractString} = nothing,
    cpu_marker::Bool = false,
)
    h = short_hash(fom_signature(cfg; tag))
    suffix = cpu_marker ? "_cpu" : ""
    return "fom_n$(cfg.n)_cfl$(fmt_token(cfg.cfl))_st$(fmt_token(cfg.stretch))_ns$(cfg.nsnapshot)_nu$(fmt_token(cfg.viscosity))$(suffix)_$(h)"
end

function rom_id(parent_fom_id::String, cfg::ROMConfig)
    h = short_hash(rom_signature(parent_fom_id, cfg))
    method = rom_method_token(cfg)
    parent_token = parent_fom_hash_token(parent_fom_id)
    basis_mode = lowercase(cfg.basis_mode)
    basis_token = if basis_mode == "global"
        "bglob"
    elseif basis_mode == "space_local"
        strat = lowercase(cfg.basis_local_strategy)
        part = lowercase(cfg.basis_partition_method)
        eta = fmt_token(cfg.basis_local_energy_fraction)
        if part in ("uniform_xy", "manual_xy")
            "bsl$(strat)_$(part)_nx$(cfg.basis_partition_nx)_ny$(cfg.basis_partition_ny)_eta$(eta)"
        else
            "bsl$(strat)_$(part)_nx$(cfg.basis_partition_nx)_eta$(eta)"
        end
    else
        "bbasis$(basis_mode)"
    end
    rank_sel = lowercase(strip(cfg.rank_selection))
    rank_token = if rank_sel in ("energy_fraction", "energyfraction", "energy", "fraction", "eta")
        "e$(fmt_token(cfg.energy_fraction))"
    else
        "r$(cfg.npod)"
    end
    lam_token = method == "trajopinf" ? "_lam$(fmt_token(cfg.lambda))" : ""
    run_tag_token = isnothing(cfg.run_tag) ? "" : "_tag$(tag_token(cfg.run_tag))"
    return "rom_from_$(parent_token)_$(method)_$(basis_token)_$(rank_token)$(lam_token)$(run_tag_token)_$(h)"
end
