using Symbolics: scalarize
using Symbolics
using ReferenceTests

@testset "stencil basics" begin

    @variables x[1:4, 1:4]
    @makearray y[1:6, 1:6] begin
        y[2:end-1, 2:end-1] => @arrayop (i, j) x[j, i]
    end
    @test isequal(scalarize(y[2,3]), x[2,1])
    @test isequal(scalarize(y[3,2]), x[1,2])
    @test_throws UndefRefError scalarize(y[1,1])

    _repr(x) = repr(SymbolicUtils.Code.toexpr(SymbolicUtils.Code.LiteralExpr(x)) |> Base.remove_linenums!)
    function test_funcs(name, f, args...)
        outplace, inplace = build_function(f, args...)

        @test_reference "build_function_tests/$name-outplace.jl" _repr(outplace)
        @test_reference "build_function_tests/$name-inplace.jl" _repr(inplace)
    end
    test_funcs("stencil-transpose-arrayop", y, x)
    # arrayop
    # constant (broadcast)
    # constant (not broadcast)
    # high-level operation

end
