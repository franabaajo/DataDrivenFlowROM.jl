if !isdefined(@__MODULE__, :_WORKFLOW_BACKEND_LOADED)
const _WORKFLOW_BACKEND_LOADED = true

using CUDA
using KernelAbstractions
import IncompressibleNavierStokes as NS

# For GPU runs, try to load CUDSS early so the extension can activate.
if CUDA.functional()
    try
        @eval using CUDSS
    catch
        # Defer hard failure to ensure_workflow_cuda_psolver!.
    end
end

function workflow_getbackend(; env_var::AbstractString = "NS_BACKEND")
    forced = lowercase(get(ENV, env_var, ""))
    if forced in ("cpu", "cpu-only", "cpu_only")
        return KernelAbstractions.CPU()
    elseif forced in ("gpu", "cuda")
        return CUDABackend()
    end

    return CUDA.functional() ? CUDABackend() : KernelAbstractions.CPU()
end

function ensure_workflow_cuda_psolver!(backend)
    if backend isa CUDA.CUDABackend
        if !CUDA.functional()
            error(
                "CUDA backend selected, but the CUDA driver is not functional / accessible on this system. " *
                "If you expected a GPU run, verify `nvidia-smi` works and that `/dev/nvidia*` exists. " *
                "To force CPU, set NS_BACKEND=cpu (or the case-specific backend env var)."
            )
        end
        try
            @eval using CUDSS
            if Base.get_extension(NS, :CUDSSExt) === nothing
                error("CUDSS extension for IncompressibleNavierStokes did not load.")
            end
        catch
            error(
                "CUDA backend selected, but CUDSS is not available. " *
                "Install and load it in Julia:\n" *
                "using Pkg; Pkg.add(\"CUDSS\"); using CUDSS"
            )
        end
    end
    return nothing
end

function workflow_backend_summary(; env_var::AbstractString = "NS_BACKEND")
    backend = workflow_getbackend(; env_var)
    if backend isa CUDA.CUDABackend
        dev = CUDA.device()
        return (backend = "gpu", device = CUDA.name(dev))
    end
    return (backend = "cpu", device = "CPU")
end

end
