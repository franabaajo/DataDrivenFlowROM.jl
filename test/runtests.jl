using TestItemRunner

@run_package_tests filter = test -> occursin(@__DIR__, test.filename)

