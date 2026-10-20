using Dates

registry_path() = joinpath(output_root(), "registry.csv")

function ensure_registry!()
    path = registry_path()
    if !isfile(path)
        open(path, "w") do io
            println(io, "kind,run_id,parent_fom_id,status,timestamp,path,n,cfl,stretch,nsnapshot,viscosity,npod,lambda,use_exact_opinf,dt_train")
        end
    end
    return path
end

function append_registry_line!(line::AbstractString)
    ensure_registry!()
    open(registry_path(), "a") do io
        println(io, line)
    end
end

function register_fom!(run_id::String, cfg::FOMConfig; status::String = "completed")
    line = join([
        "fom",
        run_id,
        "",
        status,
        string(Dates.now()),
        relpath(fom_dir(run_id), examples_root()),
        string(cfg.n),
        string(cfg.cfl),
        string(cfg.stretch),
        string(cfg.nsnapshot),
        string(cfg.viscosity),
        "",
        "",
        "",
        "",
    ], ",")
    append_registry_line!(line)
end

function register_rom!(
    run_id::String,
    parent_fom_id::String,
    romcfg::ROMConfig,
    fomcfg::FOMConfig;
    status::String = "completed",
    dt_train::Union{Nothing,Real} = nothing,
)
    use_exact_effective = romcfg.use_exact_opinf && !romcfg.use_intrusive_rom
    dt_train_effective = use_exact_effective && !isnothing(dt_train) ? string(dt_train) : ""
    line = join([
        "rom",
        run_id,
        parent_fom_id,
        status,
        string(Dates.now()),
        relpath(rom_dir(run_id), examples_root()),
        string(fomcfg.n),
        string(fomcfg.cfl),
        string(fomcfg.stretch),
        string(fomcfg.nsnapshot),
        string(fomcfg.viscosity),
        string(romcfg.npod),
        string(romcfg.lambda),
        string(use_exact_effective),
        dt_train_effective,
    ], ",")
    append_registry_line!(line)
end

function find_existing_fom(cfg::FOMConfig; cpu_marker::Bool = false, tag::Union{Nothing,AbstractString} = nothing)
    run_id = fom_id(cfg; cpu_marker = cpu_marker, tag = tag)
    return isdir(fom_dir(run_id)) ? run_id : nothing
end
find_existing_rom(parent_fom_id::String, cfg::ROMConfig) = isdir(rom_dir(rom_id(parent_fom_id, cfg))) ? rom_id(parent_fom_id, cfg) : nothing
