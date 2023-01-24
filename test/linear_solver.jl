using Symbolics
using LinearAlgebra
using Test

@variables t p(t) x y(t)
D = Differential(t)
expr = x * p + y * (1 - p) ~ 0
sol = Symbolics.solve_for(expr, p)
@test sol isa Num
@test isequal(sol, y / (y - x))
a, b, islinear = Symbolics.linear_expansion(expr, p)
@test eltype((a, b)) <: Num
@test isequal((a, b, islinear), (-(x - y), -y, true))
@test isequal(Symbolics.solve_for(x * p ~ 0, p), 0)
@test_throws Any Symbolics.solve_for(1 / x + p * p / x ~ 0, p)
@test isequal(Symbolics.solve_for(x * y ~ p, x), p / y)
@test isequal(Symbolics.solve_for(x * -y ~ p, y), -p / x)
@test isequal(Symbolics.solve_for(x * y ~ p, p), x * y)
@test isequal(Symbolics.solve_for(x^2 * y ~ p, y), p / x^2)
@test isequal(Symbolics.solve_for(x^2 * y ~ p, p), x^2 * y)
@test_throws Any Symbolics.solve_for(x^2 * y - sin(p) ~ p, p)
@test Symbolics.solve_for(x^2 * y - sin(p) ~ p, p, check = false) === nothing
@test_throws Any Symbolics.solve_for(t * D(x) ~ y, t)
@test isequal(Symbolics.solve_for(t * D(x) ~ y, D(x)), y / t)
@test isequal(Symbolics.solve_for([t * D(x) ~ y], [D(x)]), [y / t])
@test isequal(Symbolics.solve_for(t * D(x) ~ y, y), t * D(x))
Dx = D(x)
expr = Dx * x + Dx * t - 2 // 3 * x + y * Dx
a, b, islinear = Symbolics.linear_expansion(expr, x)
@test iszero(expand(a * x + b - expr))
@test isequal(Symbolics.solve_for(expr ~ Dx, Dx), (-2 // 3 * x) / (1 - t - x - y))
@test isequal(Symbolics.solve_for(expr ~ Dx, x), (t * Dx + y * Dx - Dx) / ((2 // 3) - Dx))

exprs = [3 // 2 * x + 2y + 10
         7x + 3y - 8]
xs = [x, y]
A, b, islinear = Symbolics.linear_expansion(exprs, xs)
@test islinear
@test isequal(A, [3//2 2; 7 3])
@test isequal(b, [10; -8])

@variables x y z
eqs = [2 // 1 * x + y - z ~ 2 // 1
       2 // 1 + y - z ~ 3 // 1 * x
       2 // 1 + y - 2z ~ 3 // 1 * z]
@test [2 1 -1; -3 1 -1; 0 1 -5] * Symbolics.solve_for(eqs, [x, y, z]) == [2; -2; -2]
@test isequal(Symbolics.solve_for(2 // 1 * x + y - 2 // 1 * z ~ 9 // 1 * x, 1 // 1 * x),
              (1 // 7) * (y - 2 // 1 * z))
@test isequal(Symbolics.solve_for(x + y ~ 0, x), Symbolics.solve_for([x + y ~ 0], x))
@test isequal(Symbolics.solve_for([x + y ~ 0], [x]), Symbolics.solve_for(x + y ~ 0, [x]))
@test isequal(Symbolics.solve_for(2x / z + sin(z), x), sin(z) / (-2 / z))
