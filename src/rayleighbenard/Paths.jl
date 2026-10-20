using SHA

repo_root() = normpath(joinpath(@__DIR__, "..", ".."))
examples_root() = repo_root()
rayleigh_case_name() = "RayleighBenardConvection"
default_output_root() = joinpath(repo_root(), "output", rayleigh_case_name(), "rayleighpod")
legacy_output_root() = joinpath(repo_root(), "output", "rayleighpod")

function output_root()
    env_root = get(ENV, "RAYLEIGH_OUTPUT_ROOT", "")
    if !isempty(env_root)
        return mkpath(normpath(env_root))
    end

    preferred = default_output_root()
    if isdir(preferred)
        return mkpath(preferred)
    end

    legacy = legacy_output_root()
    if isdir(legacy)
        @info "Using legacy Rayleigh-Benard output root" legacy
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
        "conductivity=$(cfg.conductivity)",
        "gravity=$(cfg.gravity)",
        "gdir=$(cfg.gdir)",
        "dodissipation=$(cfg.dodissipation)",
        "nx=$(cfg.nx)",
        "ny=$(cfg.ny)",
        "xlims=$(cfg.xlims)",
        "ylims=$(cfg.ylims)",
        "stretch_x=$(cfg.stretch_x)",
        "stretch_y=$(cfg.stretch_y)",
        "temp_hot=$(cfg.temp_hot)",
        "temp_cold=$(cfg.temp_cold)",
        "initial_temperature_profile=$(cfg.initial_temperature_profile)",
        "perturbation_amplitude=$(cfg.perturbation_amplitude)",
        "perturbation_wavenumber=$(cfg.perturbation_wavenumber)",
        "tsim=$(cfg.tsim)",
        "cfl=$(cfg.cfl)",
        "nsnapshot=$(cfg.nsnapshot)",
    ]
    if !isnothing(tag)
        push!(parts, "tag=$(tag)")
    end
    return join(parts, "|")
end

function rom_signature(parent_fom_id::String, cfg::ROMConfig)
    parts = String[
        "parent=$(parent_fom_id)",
        "npod=$(cfg.npod)",
        "use_intrusive_rom=$(cfg.use_intrusive_rom)",
        "use_exact_opinf=$(cfg.use_exact_opinf)",
        "force_explicit_euler=$(cfg.exact_opinf_force_explicit_euler)",
        "fallback_to_standard_opinf=$(cfg.fallback_to_standard_opinf)",
        "include_constant_term=$(cfg.include_constant_term)",
        "safety=$(cfg.safety)",
        "growthfac=$(cfg.growthfac)",
        "energy_fraction=$(cfg.energy_fraction)",
    ]
    if !isnothing(cfg.run_tag)
        push!(parts, "run_tag=$(cfg.run_tag)")
    end
    if cfg.training_start_time != 0.0
        push!(parts, "training_start_time=$(cfg.training_start_time)")
    end
    if !isnothing(cfg.training_end_time)
        push!(parts, "training_end_time=$(cfg.training_end_time)")
    end
    # Lambda is only used for trajectory-based operator inference.
    method = rom_method_token(cfg)
    if method == "trajopinf"
        push!(parts, "lambda=$(cfg.lambda)")
    end
    rank_sel = lowercase(strip(cfg.rank_selection))
    if !(rank_sel in ("npod", "fixed", "rank"))
        push!(parts, "rank_selection=$(rank_sel)")
    end
    if cfg.basis_centered
        push!(parts, "basis_centered=true")
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
    return join(parts, "|")
end

function switched_rom_signature(parent_fom_id::String, cfg::SwitchedROMConfig)
    validate_switched_rom_config(cfg)
    parts = String[
        "parent=$(parent_fom_id)",
        "kind=switched_rom",
        "local_rom_ids=$(join(cfg.local_rom_ids, ';'))",
        "regime_names=$(join(cfg.regime_names, ';'))",
        "start_time=$(cfg.start_time)",
        "end_time=$(cfg.end_time)",
        "initial_regime=$(cfg.initial_regime)",
        "switching_policy=$(cfg.switching_policy)",
        "min_switch_time=$(cfg.min_switch_time)",
        "pre_to_transition_threshold=$(cfg.pre_to_transition_threshold)",
        "transition_to_post_threshold=$(cfg.transition_to_post_threshold)",
        "transition_to_pre_threshold=$(cfg.transition_to_pre_threshold)",
        "post_to_transition_threshold=$(cfg.post_to_transition_threshold)",
        "allow_backward_switches=$(cfg.allow_backward_switches)",
        "pre_to_transition_time=$(cfg.pre_to_transition_time)",
        "transition_to_post_time=$(cfg.transition_to_post_time)",
        "allow_parent_fom_mismatch=$(cfg.allow_parent_fom_mismatch)",
        "safety=$(cfg.safety)",
        "growthfac=$(cfg.growthfac)",
    ]
    isnothing(cfg.run_tag) || push!(parts, "run_tag=$(cfg.run_tag)")
    return join(parts, "|")
end

function rom_method_token(cfg::ROMConfig)
    if cfg.use_intrusive_rom
        return "intrusive"
    end
    if cfg.use_exact_opinf
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
    return "fom_nx$(cfg.nx)_ny$(cfg.ny)_cfl$(fmt_token(cfg.cfl))_stx$(fmt_token(cfg.stretch_x))_sty$(fmt_token(cfg.stretch_y))_nu$(fmt_token(cfg.viscosity))_kappa$(fmt_token(cfg.conductivity))$(suffix)_$(h)"
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
    tstart_token = cfg.training_start_time == 0.0 ? "" : "_ts$(fmt_token(cfg.training_start_time))"
    tend_token = isnothing(cfg.training_end_time) ? "" : "_te$(fmt_token(cfg.training_end_time))"
    tag_token = isnothing(cfg.run_tag) ? "" : "_tag$(fmt_token(cfg.run_tag))"
    return "rom_from_$(parent_token)_$(method)_$(basis_token)_$(rank_token)$(lam_token)$(tstart_token)$(tend_token)$(tag_token)_$(h)"
end

function switched_rom_id(parent_fom_id::String, cfg::SwitchedROMConfig)
    h = short_hash(switched_rom_signature(parent_fom_id, cfg))
    parent_token = parent_fom_hash_token(parent_fom_id)
    tstart_token = cfg.start_time == 0.0 ? "" : "_ts$(fmt_token(cfg.start_time))"
    tend_token = isnothing(cfg.end_time) ? "" : "_te$(fmt_token(cfg.end_time))"
    tag_token = isnothing(cfg.run_tag) ? "" : "_tag$(fmt_token(cfg.run_tag))"
    return "rom_from_$(parent_token)_switched_lz$(tstart_token)$(tend_token)$(tag_token)_$(h)"
end
