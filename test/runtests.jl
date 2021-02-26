using SafeTestsets, Test

@safetestset "Differentiation Test" begin include("diff.jl") end
@safetestset "Array Test" begin include("array.jl") end
