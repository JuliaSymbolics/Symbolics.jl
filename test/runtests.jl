using SafeTestsets, Test

@safetestset "Differentiation Test" begin include("diff.jl") end
@safetestset "Overloading Test" begin include("overloads.jl") end
@safetestset "Build Function Test" begin include("build_function.jl") end
@safetestset "Build Function Array Test" begin include("build_function_arrayofarray.jl") end
@safetestset "Build Targets Test" begin include("build_targets.jl") end
