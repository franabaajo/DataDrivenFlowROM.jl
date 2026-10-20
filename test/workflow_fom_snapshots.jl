@testitem "Streamed snapshots validate file sizes" begin
    using JLD2

    include(joinpath(@__DIR__, "..", "src", "workflows", "FOM.jl"))

    mktempdir() do dir
        paths = snapshot_stream_paths(dir)

        len = 4
        nsnapshot = 1
        K = nsnapshot + 1

        meta = (;
            len = len,
            nsnapshot = nsnapshot,
            eltype = "Float64",
            layout = "column_major",
        )
        save_object(paths.meta, meta; iotype = IOStream)

        open(paths.u, "w") do io
            write(io, zeros(Float64, len * K))
        end
        open(paths.dudt, "w") do io
            # Write a truncated file to trigger the size check instead of SIGBUS.
            write(io, zeros(Float64, len * K - 1))
        end

        @test_throws "Snapshot stream file size mismatch" load_streamed_snapshots(dir)

        open(paths.dudt, "w") do io
            write(io, zeros(Float64, len * K))
        end

        snaps = load_streamed_snapshots(dir)
        @test size(snaps.u) == (len, K)
        @test size(snaps.dudt) == (len, K)
        @test snaps.u[1, 1] == 0.0
        @test snaps.dudt[end, end] == 0.0
    end
end

