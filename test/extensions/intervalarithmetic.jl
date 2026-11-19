using Symbolics
using IntervalArithmetic
using Test
using SymbolicUtils: operation, arguments, unwrap_const, promote_symtype

@testset "IntervalArithmetic promotion (issue #1157)" begin
    ext = Base.get_extension(Symbolics, :SymbolicsIntervalArithmeticExt)
    @test !isnothing(ext)
    @test isempty(Test.detect_ambiguities(ext))

    @test promote_type(Interval{Float64}, Num) === Num
    @test promote_type(Num, Interval{Float64}) === Num

    @variables x
    iv = interval(1, 2)

    prod1 = iv * x
    @test prod1 isa Num
    @test operation(Symbolics.value(prod1)) === *
    @test isequal(unwrap_const(arguments(Symbolics.value(prod1))[1]), iv)
    @test isequal(arguments(Symbolics.value(prod1))[2], Symbolics.value(x))

    prod2 = x * iv
    @test prod2 isa Num
    @test operation(Symbolics.value(prod2)) === *

    @test (iv + x) isa Num
    @test (x + iv) isa Num
    @test (iv - x) isa Num
    @test (x - iv) isa Num

    quot1 = iv / x
    @test quot1 isa Num
    @test operation(Symbolics.value(quot1)) === /
    @test isequal(unwrap_const(arguments(Symbolics.value(quot1))[1]), iv)
    @test isequal(arguments(Symbolics.value(quot1))[2], Symbolics.value(x))

    quot2 = x / iv
    @test quot2 isa Num
    @test operation(Symbolics.value(quot2)) === /

    # Backslash is Div with swapped arguments.
    bslash1 = iv \ x
    @test bslash1 isa Num
    @test operation(Symbolics.value(bslash1)) === /
    @test isequal(arguments(Symbolics.value(bslash1))[1], Symbolics.value(x))
    @test isequal(unwrap_const(arguments(Symbolics.value(bslash1))[2]), iv)

    bslash2 = x \ iv
    @test bslash2 isa Num
    @test operation(Symbolics.value(bslash2)) === /

    # Thin unit interval should cancel under multiplication.
    @test isequal(interval(1, 1) * x, x)

    # Display must not throw (SymbolicUtils `_isunit` / `iszero` use `==` by default).
    @test sprint(show, prod1) isa String
    @test sprint(show, quot1) isa String
    @test sprint(show, quot2) isa String
    @test sprint(show, interval(-1, 1) * x) isa String

    # Interval×Interval substitute must not hit promote_symtype self-ambiguity.
    subbed = substitute(x / iv, Dict(x => interval(4, 8)))
    @test subbed isa Num
    @test isequal(subbed, interval(2, 8))

    # Complex symtypes must fall through to SymbolicUtils' Complex promote_symtype.
    @test promote_symtype(/, Complex{Real}, Interval{Float64}) === Complex{Real}
    @test promote_symtype(/, Interval{Float64}, Complex{Real}) === Complex{Real}

    # simplify on interval-coefficient quotients is intentionally unsupported
    # (would need further Base piracy: isinteger, Int/Integer/Float64 conversion).
    # `(iv*x)^interval(2,2)` constructs, but printing / Float64 conversion throws —
    # documented in the PR body rather than asserted here.
    @test_throws IntervalArithmetic.InconclusiveBooleanOperation simplify(2x / iv)
    @test_throws IntervalArithmetic.InconclusiveBooleanOperation simplify(iv * x / 2)
    @test_throws MethodError simplify(2x / interval(2, 2))
end
