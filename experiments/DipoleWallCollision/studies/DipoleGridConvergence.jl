# # Two-dimensional dipole grid convergence
#
# This script evaluates how the full-order model changes with grid resolution for
# the dipole setup. It compares the total kinetic energy and enstrophy histories
# for a set of grids against a finest-grid reference solution.

using CairoMakie
using CUDA
using CUDSS
using JLD2
using ProgressMeter

import IncompressibleNavierStokes as NS

const INTERPOLATE_OMEGA_P = getproperty(NS, Symbol("interpolate_\u03c9_p"))
const SETUP_D = Symbol("\u0394")
const SETUP_DU = Symbol("\u0394u")

# Define all parameters here for convenience
params() = (;
    ## Fluid
    viscosity = 1 / 2500, # Viscosity

    ## Initial conditions
    r0 = 0.1,             # Initial vortex radius
    p1 = (0.0, 0.1),      # Position of first vortex
    p2 = (0.0, -0.1),     # Position of second vortex
    initialenergy = 2.0,  # Initial kinetic energy

    ## Grid
    stretch = 1.0,        # Grid stretching factor (`nothing` for uniform grid)
    lims = (-1.0, 1.0),   # Domain limits
    nlist = [128, 256, 512, 1024],
    nref = 1024,

    ## Full-order model
    tsim = 1.0,           # Simulation time
    cfl = 0.9,            # CFL number
    nsnapshot = 100,      # How many snapshots to save (in addition to the initial snapshot)
)

# Output directory (for saving plots and diagnostics)
output() = mkpath(joinpath(@__DIR__, "..", "output", "DipoleGridConvergence"))

# If a CUDA-compatible GPU is available, use it. Otherwise, use CPU.
getbackend() =
    if CUDA.functional()
        CUDABackend()
    else
        NS.KernelAbstractions.CPU()
    end

# Domain setup
function getsetup(n)
    (; lims, stretch) = params()
    ax = if isnothing(stretch)
        range(lims..., n + 1)
    else
        NS.tanh_grid(lims..., n, stretch)
    end
    NS.Setup(;
        x = (ax, ax),
        boundary_conditions = (;
            u = (
                (NS.DirichletBC(), NS.DirichletBC()),
                (NS.DirichletBC(), NS.DirichletBC()),
            ),
        ),
        backend = getbackend(),
    )
end

# Initial conditions: Two vortices (dipole)
function U(dim, x, y)
    @inline
    (; p1, p2, r0) = params()
    x1, y1 = p1
    x2, y2 = p2
    r1 = (x - x1)^2 + (y - y1)^2
    r2 = (x - x2)^2 + (y - y2)^2
    d1 = ifelse(dim == 1, -(y - y1), x - x1)
    d2 = ifelse(dim == 1, y - y2, -(x - x2))
    return d1 * exp(-r1 / r0^2) + d2 * exp(-r2 / r0^2)
end

# Scale velocity to have the correct initial kinetic energy
function scalevelocity!(u, setup)
    (; initialenergy) = params()
    u1 = selectdim(u, 3, 1)
    u2 = selectdim(u, 3, 2)
    d = getproperty(setup, SETUP_D)
    du = getproperty(setup, SETUP_DU)

    d1x = du[1]
    d1y = d[2]'
    d2x = d[1]
    d2y = du[2]'

    eu = @. d1x * d1y * u1^2 / 2
    ev = @. d2x * d2y * u2^2 / 2
    kin = sum(view(eu, setup.Iu[1])) + sum(view(ev, setup.Iu[2]))

    @. u = sqrt(initialenergy / kin) * u
    return
end

function total_enstrophy(u, setup)
    omega = NS.vorticity(u, setup)
    omega_p = INTERPOLATE_OMEGA_P(omega, setup)
    enstrophy_density = @. omega_p^2 / 2
    enstrophy_density = NS.scalewithvolume(enstrophy_density, setup)
    sum(view(enstrophy_density, setup.Ip))
end

function relative_error(history, reference; floor = 1e-12)
    @assert length(history) == length(reference)
    denom = max(maximum(abs, reference), floor)
    abs.(history .- reference) ./ denom
end

function run_diagnostics(n)
    (; viscosity, tsim, cfl, nsnapshot) = params()

    setup = getsetup(n)
    psolver = NS.default_psolver(setup)
    u = NS.velocityfield(setup, U; psolver)
    scalevelocity!(u, setup)

    state = (; u = copy(u))
    method = NS.LMWray3()
    ode_cache = NS.get_cache(method, state, setup)
    force_cache = NS.get_cache(NS.navierstokes!, setup)
    stepper = NS.create_stepper(method; setup, psolver, state, t = 0.0)

    times = zeros(nsnapshot + 1)
    kinetic_energy = zeros(nsnapshot + 1)
    enstrophy = zeros(nsnapshot + 1)

    prog = Progress(nsnapshot + 1; desc = "Running dipole diagnostics on n = $n")
    for isnap = 0:nsnapshot
        ## Step until next save point.
        ## For the first step, the while loop is never entered.
        tstop = isnap / nsnapshot * tsim
        isubstep = 0
        while stepper.t < prevfloat(tstop)
            ## Change timestep based on operators
            Δt =
                cfl *
                NS.propose_timestep(NS.navierstokes!, stepper.state, setup, (; viscosity))

            ## Make sure not to step past `tstop`
            Δt = min(Δt, tstop - stepper.t)

            ## Perform a single time step with the time integration method
            stepper = NS.timestep!(
                method,
                NS.navierstokes!,
                stepper,
                Δt;
                params = (; viscosity),
                ode_cache,
                force_cache,
            )
            isubstep += 1
        end

        times[isnap + 1] = stepper.t
        kinetic_energy[isnap + 1] = NS.total_kinetic_energy(stepper.state.u, setup)
        enstrophy[isnap + 1] = total_enstrophy(stepper.state.u, setup)

        next!(prog; showvalues = [("substeps", isubstep), ("time", stepper.t)])
    end

    return (; n, times, kinetic_energy, enstrophy)
end

function summarize_errors(results, reference)
    summaries = NamedTuple[]
    for result in results
        if result.n == reference.n
            continue
        end
        ke_error = relative_error(result.kinetic_energy, reference.kinetic_energy)
        en_error = relative_error(result.enstrophy, reference.enstrophy)
        push!(
            summaries,
            (
                n = result.n,
                max_ke_rel_error = maximum(ke_error),
                final_ke_rel_error = ke_error[end],
                max_en_rel_error = maximum(en_error),
                final_en_rel_error = en_error[end],
            ),
        )
    end
    summaries
end

function save_plots(results, reference)
    fig = Figure(size = (1000, 900))

    ax_ke = Axis(
        fig[1, 1];
        title = "Total kinetic energy",
        xlabel = "t",
        ylabel = "E(t)",
    )
    ax_en = Axis(
        fig[2, 1];
        title = "Total enstrophy",
        xlabel = "t",
        ylabel = "Z(t)",
    )
    ax_ke_err = Axis(
        fig[3, 1];
        title = "Relative kinetic-energy error vs n = $(reference.n)",
        xlabel = "t",
        ylabel = "relative error",
        yscale = log10,
    )
    ax_en_err = Axis(
        fig[4, 1];
        title = "Relative enstrophy error vs n = $(reference.n)",
        xlabel = "t",
        ylabel = "relative error",
        yscale = log10,
    )

    colors = Makie.wong_colors()
    for (i, result) in enumerate(results)
        color = colors[mod1(i, length(colors))]
        label = "n = $(result.n)"
        lines!(ax_ke, result.times, result.kinetic_energy; label, color)
        lines!(ax_en, result.times, result.enstrophy; label, color)

        if result.n != reference.n
            ke_error = relative_error(result.kinetic_energy, reference.kinetic_energy)
            en_error = relative_error(result.enstrophy, reference.enstrophy)
            lines!(ax_ke_err, result.times, ke_error; label, color)
            lines!(ax_en_err, result.times, en_error; label, color)
        end
    end

    axislegend(ax_ke; position = :rb)
    axislegend(ax_en; position = :rb)
    axislegend(ax_ke_err; position = :rb)
    axislegend(ax_en_err; position = :rb)

    save("$(output())/diagnostics.png", fig)
    return fig
end

function save_summary(results, summaries)
    save_object("$(output())/results.jld2", (; params = params(), results, summaries))

    open("$(output())/summary.txt", "w") do io
        println(io, "Dipole grid convergence study")
        println(io, "Reference grid: n = $(params().nref)")
        println(io, "Snapshot count: $(params().nsnapshot + 1)")
        println(io)
        println(io, "Relative errors are computed against the n = $(params().nref) history.")
        println(io)
        for summary in summaries
            println(io, "n = $(summary.n)")
            println(io, "  max kinetic-energy relative error = $(summary.max_ke_rel_error)")
            println(io, "  final kinetic-energy relative error = $(summary.final_ke_rel_error)")
            println(io, "  max enstrophy relative error = $(summary.max_en_rel_error)")
            println(io, "  final enstrophy relative error = $(summary.final_en_rel_error)")
            println(io)
        end
    end
end

function main()
    (; nlist, nref) = params()
    @assert nref in nlist "The reference grid must be included in nlist."

    results = [run_diagnostics(n) for n in nlist]
    reference = only(filter(result -> result.n == nref, results))
    summaries = summarize_errors(results, reference)

    save_summary(results, summaries)
    save_plots(results, reference)

    for summary in summaries
        @info(
            "Grid comparison against reference",
            n = summary.n,
            max_ke_rel_error = summary.max_ke_rel_error,
            final_ke_rel_error = summary.final_ke_rel_error,
            max_en_rel_error = summary.max_en_rel_error,
            final_en_rel_error = summary.final_en_rel_error,
        )
    end
end

main()
