using SafeTestsets, Test, Pkg
using Symbolics

const GROUP = get(ENV, "GROUP", "All")

function activate_downstream_env()
    Pkg.activate("downstream")
    Pkg.develop(PackageSpec(path=dirname(@__DIR__)))
    Pkg.instantiate()
end

if haskey(ENV, "BENCHMARK_ONLY")
    include("benchmark.jl")
end

# this needs to be defined at top level
limit(a, N) = a == N + 1 ? 1 : a == 0 ? N : a
@register_symbolic limit(a, N)::Integer

if GROUP == "All" || GROUP == "Core"
    @safetestset "Quality Assurance" begin include("qa.jl") end
    @safetestset "Macro Test" begin include("macro.jl") end
    @safetestset "Arrays" begin include("arrays.jl") end
    @safetestset "View-setting" begin include("stencils.jl") end
    @safetestset "Complex" begin include("complex.jl") end
    @safetestset "Semi-polynomial" begin include("semipoly.jl") end
    @safetestset "Fuzz Arrays" begin include("fuzz-arrays.jl") end
    @safetestset "Differentiation Test" begin include("diff.jl") end
    @safetestset "Difference Test" begin include("difference.jl") end
    @safetestset "Degree Test" begin include("degree.jl") end
    @safetestset "Coeff Test" begin include("coeff.jl") end
    @safetestset "Parsing Test" begin include("parsing.jl") end
    @safetestset "Is Linear or Affine Test" begin include("islinear_affine.jl") end
    @safetestset "Linear Solver Test" begin include("linear_solver.jl") end
    @safetestset "Algebraic Solver Test" begin include("solver.jl") end
    @safetestset "Groebner Bases Test" begin include("groebner_basis.jl") end
    @safetestset "Overloading Test" begin include("overloads.jl") end
    @safetestset "Build Function Test" begin include("build_function.jl") end
    @safetestset "Build Function Array Test" begin include("build_function_arrayofarray.jl") end
    @safetestset "Build Function Array Test Named Tuples" begin include("build_function_arrayofarray_named_tuples.jl") end
    VERSION >= v"1.9" && @safetestset "Build Targets Test" begin include("build_targets.jl") end
    @safetestset "Latexify Test" begin include("latexify.jl") end
    @safetestset "Domain Test" begin include("domains.jl") end
    @safetestset "SymPy Test" begin include("sympy.jl") end
    @safetestset "Inequality Test" begin include("inequality.jl") end
    @safetestset "Integral Test" begin include("integral.jl") end
    @safetestset "CartesianIndex Test" begin include("cartesianindex.jl") end
    @safetestset "LogExpFunctions Test" begin include("logexpfunctions.jl") end
end

if GROUP == "All" || GROUP == "Core" || GROUP == "SymbolicIndexingInterface"
    @safetestset "SymbolicIndexingInterface Trait Test" begin
        include("symbolic_indexing_interface_trait.jl")
    end
    @safetestset "SymbolicIndexingInterface Parameter Indexing Test" begin
        include("symbolic_indexing_interface_parameter_indexing.jl")
    end
end

if GROUP == "Downstream"
    activate_downstream_env()
    #@time @safetestset "ParameterizedFunctions MATLABDiffEq Regression Test" begin include("downstream/ParameterizedFunctions_MATLAB.jl") end
end
