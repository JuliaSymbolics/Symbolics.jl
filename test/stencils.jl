using Symbolics: scalarize
using Symbolics
using ReferenceTests
using Test
using SymbolicUtils.Code: toexpr, LiteralExpr

_repr(x) = repr(toexpr(LiteralExpr(x)) |> Base.remove_linenums!)
function test_funcs(name, f, args...; broken=false)
    outplace, inplace = build_function(f, args...)
    if broken
        @test_broken open(x->read(x, String), "build_function_tests/$name-outplace.jl") == _repr(outplace)
        @test_broken open(x->read(x, String), "build_function_tests/$name-inplace.jl") == _repr(outplace)
    else
        @test_reference "build_function_tests/$name-outplace.jl" _repr(outplace)
        @test_reference "build_function_tests/$name-inplace.jl" _repr(inplace)
    end
end

@testset "array codegen basics" begin
    @variables x[1:4, 1:4]

    broken = VERSION < v"1.7"
    # Simple test with ArrayOp and no term
    test_funcs("transpose", @arrayop((i, j), x[j, i]), x)

    # Simple test with ArrayOp and Term
    test_funcs("transpose-term", @arrayop((i, j), x[j, i], term=x'), x)

    # Partial view set to an arrayop
    @makearray y[1:6, 1:6] begin
        y[2:end-1, 2:end-1] => @arrayop (i, j) x[j, i]
    end
    @test isequal(scalarize(y[2,3]), x[2,1])
    @test isequal(scalarize(y[3,2]), x[1,2])

    # Test UndefRef is thrown
    @test_throws UndefRefError scalarize(y[1,1])

    test_funcs("stencil-transpose-arrayop", y, x)

    # Fill zero
    @makearray y[1:6, 1:6] begin
        y[:, :] => 0
        y[2:end-1, 2:end-1] => x .+ x' .+ 1
    end

    @test iszero(scalarize(y[1,1]))

    test_funcs("stencil-broadcast", y, x)

    @variables x[1:5, 1:5]
    @makearray y[1:5, 1:5] begin
        y[:, :] => 0
        y[2:end-1, 2:end-1] => @arrayop (i, j) (x[i+1,j] + x[i-1, j] + x[i, j+1] + x[i, j-1])/2
    end

    @test iszero(scalarize(y[1,1]))
    test_funcs("stencil-extents", y, x, broken=broken)

    @variables u[1:5, 1:5]
    n = 5
    limit = Main.limit
    y = @arrayop (i, j) u[limit(i-1, n), limit(j+1,n)] i in 1:n j in 1:n
    test_funcs("manual-limits", y, u)

    z = @arrayop (i, j) y[j, i]
    test_funcs("intermediate-exprs", z, u)
end
