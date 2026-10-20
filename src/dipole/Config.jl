using TOML

Base.@kwdef struct FOMConfig
    viscosity::Float64 = 1 / 2500
    r0::Float64 = 0.1
    p1::NTuple{2,Float64} = (0.0, 0.1)
    p2::NTuple{2,Float64} = (0.0, -0.1)
    initialenergy::Float64 = 2.0

    n::Int = 768
    stretch::Union{Nothing,Float64} = 1.0
    lims::NTuple{2,Float64} = (-1.0, 1.0)

    tsim::Float64 = 1.0
    cfl::Float64 = 0.67
    nsnapshot::Int = 100
end

Base.@kwdef struct ROMConfig
    npod::Int = 20
    # Rank selection strategy for global POD basis:
    # - "npod": use `npod` (fixed rank)
    # - "energy_fraction": choose the smallest rank capturing `energy_fraction`
    rank_selection::String = "npod"
    lambda::Float64 = 1e-2
    safety::Float64 = 0.03
    growthfac::Float64 = 1.05
    energy_fraction::Float64 = 0.99
    use_intrusive_rom::Bool = false
    use_exact_opinf::Bool = true
    exact_opinf_force_explicit_euler::Bool = true
    fallback_to_standard_opinf::Bool = true
    run_tag::Union{Nothing,String} = nothing

    # POD basis mode:
    # - "global": standard POD (current default)
    # - "space_local": space-local POD on row blocks
    basis_mode::String = "global"
    basis_local_strategy::String = "S3" # "S1" (fixed total rank) or "S3" (per-block energy)
    basis_partition_method::String = "uniform_x" # "uniform_x", "uniform_xy", "manual_x", or "manual_xy"
    basis_partition_nx::Int = 1
    basis_partition_ny::Int = 1
    basis_partition_xcuts::Union{Nothing,Vector{Float64}} = nothing
    basis_partition_ycuts::Union{Nothing,Vector{Float64}} = nothing
    basis_local_energy_fraction::Float64 = 0.99
    basis_local_ranks::Union{Nothing,Vector{Int}} = nothing
    basis_space_local_sparse::Bool = true
end

_stretch_to_toml(x) = isnothing(x) ? "uniform" : x
_stretch_from_toml(x) = x == "uniform" ? nothing : Float64(x)

function fom_to_dict(cfg::FOMConfig)
    Dict(
        "viscosity" => cfg.viscosity,
        "r0" => cfg.r0,
        "p1" => collect(cfg.p1),
        "p2" => collect(cfg.p2),
        "initialenergy" => cfg.initialenergy,
        "n" => cfg.n,
        "stretch" => _stretch_to_toml(cfg.stretch),
        "lims" => collect(cfg.lims),
        "tsim" => cfg.tsim,
        "cfl" => cfg.cfl,
        "nsnapshot" => cfg.nsnapshot,
    )
end

function rom_to_dict(cfg::ROMConfig)
    d = Dict{String,Any}(
        "npod" => cfg.npod,
        "rank_selection" => cfg.rank_selection,
        "lambda" => cfg.lambda,
        "safety" => cfg.safety,
        "growthfac" => cfg.growthfac,
        "energy_fraction" => cfg.energy_fraction,
        "use_intrusive_rom" => cfg.use_intrusive_rom,
        "use_exact_opinf" => cfg.use_exact_opinf,
        "exact_opinf_force_explicit_euler" => cfg.exact_opinf_force_explicit_euler,
        "fallback_to_standard_opinf" => cfg.fallback_to_standard_opinf,
        "basis_mode" => cfg.basis_mode,
        "basis_local_strategy" => cfg.basis_local_strategy,
        "basis_partition_method" => cfg.basis_partition_method,
        "basis_partition_nx" => cfg.basis_partition_nx,
        "basis_partition_ny" => cfg.basis_partition_ny,
        "basis_local_energy_fraction" => cfg.basis_local_energy_fraction,
        "basis_space_local_sparse" => cfg.basis_space_local_sparse,
    )

    # TOML does not support null. Only write optional fields if present.
    if !isnothing(cfg.basis_partition_xcuts)
        d["basis_partition_xcuts"] = cfg.basis_partition_xcuts
    end
    if !isnothing(cfg.basis_partition_ycuts)
        d["basis_partition_ycuts"] = cfg.basis_partition_ycuts
    end
    if !isnothing(cfg.basis_local_ranks)
        d["basis_local_ranks"] = cfg.basis_local_ranks
    end
    if !isnothing(cfg.run_tag)
        d["run_tag"] = cfg.run_tag
    end

    return d
end

function dict_to_fom(d::AbstractDict)
    FOMConfig(
        viscosity = Float64(d["viscosity"]),
        r0 = Float64(d["r0"]),
        p1 = Tuple(Float64.(d["p1"])),
        p2 = Tuple(Float64.(d["p2"])),
        initialenergy = Float64(d["initialenergy"]),
        n = Int(d["n"]),
        stretch = _stretch_from_toml(d["stretch"]),
        lims = Tuple(Float64.(d["lims"])),
        tsim = Float64(d["tsim"]),
        cfl = Float64(d["cfl"]),
        nsnapshot = Int(d["nsnapshot"]),
    )
end

function dict_to_rom(d::AbstractDict)
    ROMConfig(
        npod = Int(d["npod"]),
        rank_selection = haskey(d, "rank_selection") ? String(d["rank_selection"]) : "npod",
        lambda = Float64(d["lambda"]),
        safety = Float64(d["safety"]),
        growthfac = Float64(d["growthfac"]),
        energy_fraction = haskey(d, "energy_fraction") ? Float64(d["energy_fraction"]) : 0.99,
        use_intrusive_rom = haskey(d, "use_intrusive_rom") ? Bool(d["use_intrusive_rom"]) : false,
        use_exact_opinf = haskey(d, "use_exact_opinf") ? Bool(d["use_exact_opinf"]) : true,
        exact_opinf_force_explicit_euler = haskey(d, "exact_opinf_force_explicit_euler") ?
            Bool(d["exact_opinf_force_explicit_euler"]) : true,
        fallback_to_standard_opinf = haskey(d, "fallback_to_standard_opinf") ?
            Bool(d["fallback_to_standard_opinf"]) : true,
        run_tag = haskey(d, "run_tag") && !isempty(strip(String(d["run_tag"]))) ?
            strip(String(d["run_tag"])) : nothing,

        basis_mode = haskey(d, "basis_mode") ? String(d["basis_mode"]) : "global",
        basis_local_strategy = haskey(d, "basis_local_strategy") ? String(d["basis_local_strategy"]) : "S3",
        basis_partition_method = haskey(d, "basis_partition_method") ? String(d["basis_partition_method"]) : "uniform_x",
        basis_partition_nx = haskey(d, "basis_partition_nx") ? Int(d["basis_partition_nx"]) : 1,
        basis_partition_ny = haskey(d, "basis_partition_ny") ? Int(d["basis_partition_ny"]) : 1,
        basis_partition_xcuts = haskey(d, "basis_partition_xcuts") && !isnothing(d["basis_partition_xcuts"]) ?
            Float64.(d["basis_partition_xcuts"]) : nothing,
        basis_partition_ycuts = haskey(d, "basis_partition_ycuts") && !isnothing(d["basis_partition_ycuts"]) ?
            Float64.(d["basis_partition_ycuts"]) : nothing,
        basis_local_energy_fraction = haskey(d, "basis_local_energy_fraction") ? Float64(d["basis_local_energy_fraction"]) : 0.99,
        basis_local_ranks = haskey(d, "basis_local_ranks") && !isnothing(d["basis_local_ranks"]) ?
            Int.(d["basis_local_ranks"]) : nothing,
        basis_space_local_sparse = haskey(d, "basis_space_local_sparse") ? Bool(d["basis_space_local_sparse"]) : true,
    )
end

function save_fom_config(path::AbstractString, cfg::FOMConfig)
    open(path, "w") do io
        TOML.print(io, fom_to_dict(cfg))
    end
end

function save_rom_config(path::AbstractString, cfg::ROMConfig)
    open(path, "w") do io
        TOML.print(io, rom_to_dict(cfg))
    end
end

load_fom_config(path::AbstractString) = dict_to_fom(TOML.parsefile(path))
load_rom_config(path::AbstractString) = dict_to_rom(TOML.parsefile(path))
