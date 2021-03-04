using SafeTestsets, Test

function activate_downstream_env()
    Pkg.activate("downstream")
    Pkg.develop(PackageSpec(path=dirname(@__DIR__)))
    Pkg.instantiate()
end

if !is_APPVEYOR && GROUP == "Core"
    @safetestset "Differentiation Test" begin include("diff.jl") end
    @safetestset "Overloading Test" begin include("overloads.jl") end
    @safetestset "Build Function Test" begin include("build_function.jl") end
    @safetestset "Build Function Array Test" begin include("build_function_arrayofarray.jl") end
    @safetestset "Build Targets Test" begin include("build_targets.jl") end
end

if !is_APPVEYOR && GROUP == "Downstream"
    activate_downstream_env()
    @time @safetestset "Unitful" begin include("downstream/ParameterizedFunctions_MATLAB.jl") end
end
