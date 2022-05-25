using Symbolics: scalarize
using Symbolics
using ReferenceTests
using Test
using SymbolicUtils.Code: toexpr, LiteralExpr

_repr(x) = repr(toexpr(LiteralExpr(x)) |> Base.remove_linenums!)
function test_funcs(name, f, args...)
    outplace, inplace = build_function(f, args...)

    @test_reference "build_function_tests/$name-outplace.jl" _repr(outplace)
    @test_reference "build_function_tests/$name-inplace.jl" _repr(inplace)
end

@testset "array codegen basics" begin
    @variables x[1:4, 1:4]
    test_funcs("transpose", @arrayop((i, j), x[j, i]), x)
    test_funcs("transpose-term", @arrayop((i, j), x[j, i], term=x'), x)
    @makearray y[1:6, 1:6] begin
        y[2:end-1, 2:end-1] => @arrayop (i, j) x[j, i]
    end
    @test isequal(scalarize(y[2,3]), x[2,1])
    @test isequal(scalarize(y[3,2]), x[1,2])
    @test_throws UndefRefError scalarize(y[1,1])

    test_funcs("stencil-transpose-arrayop", y, x)

    @makearray y[1:6, 1:6] begin
        y[:, :] => 0
        y[2:end-1, 2:end-1] => x .+ x' .+ 1
    end

    test_funcs("stencil-broadcast", y, x)

    @variables x[1:5, 1:5]
    @makearray y[1:5, 1:5] begin
        y[:, :] => 0
        y[2:end-1, 2:end-1] => @arrayop (i, j) (x[i+1,j] + x[i-1, j] + x[i, j+1] + x[i, j-1])/2
    end

    test_funcs("stencil-extents", y, x)

end
