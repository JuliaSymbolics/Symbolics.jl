using Test, BenchmarkTools

include(joinpath(@__DIR__, "..", "benchmark", "benchmarks.jl"))

@test run(SUITE["dstar_jacobian"]["simple/jacobian"]; samples = 1, evals = 1) isa BenchmarkTools.Trial
