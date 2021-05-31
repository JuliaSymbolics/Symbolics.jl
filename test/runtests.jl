using SafeTestsets, Test, Pkg

const GROUP = get(ENV, "GROUP", "All")

function activate_downstream_env()
    Pkg.activate("downstream")
    Pkg.develop(PackageSpec(path=dirname(@__DIR__)))
    Pkg.instantiate()
end

if GROUP == "All" || GROUP == "Core"
    @safetestset "Macro Test" begin include("macro.jl") end
    @safetestset "Differentiation Test" begin include("diff.jl") end
    @safetestset "Is Linear or Affine Test" begin include("islinear_affine.jl") end
    @safetestset "Linear Solver Test" begin include("linear_solver.jl") end
    @safetestset "Overloading Test" begin include("overloads.jl") end
    @safetestset "Build Function Test" begin include("build_function.jl") end
    @safetestset "Build Function Array Test" begin include("build_function_arrayofarray.jl") end
    @safetestset "Build Targets Test" begin include("build_targets.jl") end
    @safetestset "Latexify Test" begin include("latexify.jl") end
end

if GROUP == "Downstream"
    activate_downstream_env()
    #@time @safetestset "ParameterizedFunctions MATLABDiffEq Regression Test" begin include("downstream/ParameterizedFunctions_MATLAB.jl") end
end
