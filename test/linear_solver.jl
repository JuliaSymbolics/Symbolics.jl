using Symbolics
using LinearAlgebra
using Test

@variables t p(t) x y(t)
D = Differential(t)

# Test division by zero gives a NaN
# https://github.com/JuliaSymbolics/Symbolics.jl/issues/1540
@test Symbolics.solve_for((x~0), y) === NaN

expr = x * p + y * (1 - p) ~ 0
sol = Symbolics.solve_for(expr, p)
@test sol isa Num
@test isequal(sol, y/(y - x))
a, b, islinear = Symbolics.linear_expansion(expr, p)
@test eltype((a, b)) <: Num
@test isequal((a, b, islinear), (-(x - y), -y, true))
@test isequal(Symbolics.solve_for(x * p ~ 0, p), 0)
@test_throws Any Symbolics.solve_for(1/x + p * p/x ~ 0, p)
@test isequal(Symbolics.solve_for(x * y ~ p, x), p / y)
@test isequal(Symbolics.solve_for(x * -y ~ p, y), -p / x)
@test isequal(Symbolics.solve_for(x * y ~ p, p), x * y)
@test isequal(Symbolics.solve_for(x^2 * y ~ p, y), p / x^2)
@test isequal(Symbolics.solve_for(x^2 * y ~ p, p), x^2 * y)
@test_throws Any Symbolics.solve_for(x^2 * y - sin(p) ~ p, p)
@test Symbolics.solve_for(x^2 * y - sin(p) ~ p, p, check=false) === nothing
@test_throws Any Symbolics.solve_for(t*D(x) ~ y, t)
@test isequal(Symbolics.solve_for(t*D(x) ~ y, D(x)), y/t)
@test isequal(Symbolics.solve_for([t*D(x) ~ y], [D(x)]), [y/t])
@test isequal(Symbolics.solve_for(t*D(x) ~ y, y), t*D(x))
Dx = D(x)
expr = Dx * x + Dx*t - 2//3*x + y*Dx
a, b, islinear = Symbolics.linear_expansion(expr, x)
@test iszero(expand(a * x + b - expr))
@test isequal(Symbolics.solve_for(expr ~ Dx, Dx), (-2//3*x)/(1 - t - x - y))
@test isequal(Symbolics.solve_for(expr ~ Dx, x), (t*Dx + y*Dx - Dx) / ((2//3) - Dx))

exprs = [
 3//2*x + 2y + 10
 7x + 3y - 8
]
xs = [x, y]
A, b, islinear = Symbolics.linear_expansion(exprs, xs)
@test islinear
@test isequal(A, [3//2 2; 7 3])
@test isequal(b, [10; -8])

@variables x y z
eqs = [
        2//1 * x + y - z ~ 2//1
        2//1 + y - z ~ 3//1*x
        2//1 + y - 2z ~ 3//1*z
      ]
@test [2 1 -1; -3 1 -1; 0 1 -5] * Symbolics.solve_for(eqs, [x, y, z]) == [2; -2; -2]
@test isequal(Symbolics.solve_for(2//1*x + y - 2//1*z ~ 9//1*x, 1//1*x), (1//7)*(y - 2//1*z))
@test isequal(Symbolics.solve_for(x + y ~ 0, x), Symbolics.solve_for([x + y ~ 0], x))
@test isequal(Symbolics.solve_for([x + y ~ 0], [x]), Symbolics.solve_for(x + y ~ 0, [x]))
@test isequal(Symbolics.solve_for(2x/z + sin(z), x), sin(z) / (-2 / z))

@variables t x
D = Symbolics.Difference(t; dt=1)
a, b, islinear = Symbolics.linear_expansion(D(x) - x, x)
@test islinear
@test isequal(a, -1)
@test isequal(b, D(x))

@testset "linear_expansion with array variables" begin
    @variables x[1:2] y[1:2] z(..)
    @test !Symbolics.linear_expansion(z(x) + x[1], x[1])[3]
    @test !Symbolics.linear_expansion(z(x[1]) + x[1], x[1])[3]
    a, b, islin = Symbolics.linear_expansion(z(x[2]) + x[1], x[1])
    @test islin && isequal(a, 1) && isequal(b, z(x[2]))
    a, b, islin = Symbolics.linear_expansion((x + x)[1], x[1])
    @test islin && isequal(a, 2) && isequal(b, 0)
    a, b, islin = Symbolics.linear_expansion(y[1], x[1])
    @test islin && isequal(a, 0) && isequal(b, y[1])
    @test !Symbolics.linear_expansion(z([x...]), x[1])[3]
    @test !Symbolics.linear_expansion(z(collect(Symbolics.unwrap(x))), x[1])[3]
    @test !Symbolics.linear_expansion(z([x, 2x]), x[1])[3]

    @variables x[0:2]
    a, b, islin = Symbolics.linear_expansion(x[0] - z(x[1]), z(x[1]))
    @test islin && isequal(a, -1) && isequal(b, x[0])

    @variables x::Vector{Real}
    a, b, islin = Symbolics.linear_expansion(x[0] - z(x[1]), z(x[1]))
    @test islin && isequal(a, -1) && isequal(b, x[0])

    @variables x y
    @test !Symbolics.linear_expansion(x + z([x, y]), y)[3]
end

@testset "linear_expansion of ifelse" begin
    @variables p
    a, b, islin = Symbolics.linear_expansion(ifelse(p < 1, 2p + 1, 3p + 2), p)
    @test islin
    @test isequal(a, ifelse(p < 1, 2, 3))
    @test isequal(b, ifelse(p < 1, 1, 2))

    a, b, islin = Symbolics.linear_expansion(ifelse(p < 1, 2p, 3p), p)
    @test islin
    @test isequal(a, ifelse(p < 1, 2, 3))
    @test isequal(b, 0)

    a, b, islin = Symbolics.linear_expansion(ifelse(p < 1, 2p, 2p + 1), p)
    @test islin
    @test isequal(a, 2)
    @test isequal(b, ifelse(p < 1, 0, 1))

    @test !Symbolics.linear_expansion(ifelse(p < 1, p^2, 2p), p)[3]
    @test !Symbolics.linear_expansion(ifelse(p < 1, 2p, p^2), p)[3]
end
