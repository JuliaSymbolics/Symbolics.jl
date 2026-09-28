using Symbolics
using IntervalArithmetic
using Test
using SymbolicUtils: operation, arguments, unwrap_const

@testset "IntervalArithmetic promotion (issue #1157)" begin
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

    # Thin unit interval should cancel under multiplication.
    @test isequal(interval(1, 1) * x, x)

    # Display must not throw (SymbolicUtils `_isunit` / `iszero` use `==` by default).
    @test sprint(show, prod1) isa String
    @test sprint(show, quot1) isa String
    @test sprint(show, quot2) isa String
    @test sprint(show, interval(-1, 1) * x) isa String
end
