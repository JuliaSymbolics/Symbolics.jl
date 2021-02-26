using SafeTestsets, Test

@safetestset "Differentiation Test" begin include("diff.jl") end
@safetestset "Overloading Test" begin include("overloads.jl") end

if haskey(ENV, "CI")
    using Pkg
    Pkg.add(url="https://github.com/SciML/ModelingToolkit.jl.git", rev="master")
    Pkg.test("ModelingToolkit")
end
