using Symbolics
using Symbolics: value, unwrap
import SymbolicUtils as SU
using LinearAlgebra
using Test

@variables t p(t) x y(t)
D = Differential(t)

# Test division by zero gives a NaN
# https://github.com/JuliaSymbolics/Symbolics.jl/issues/1540
@test value(Symbolics.symbolic_linear_solve((x~0), y)) === NaN

expr = x * p + y * (1 - p) ~ 0
sol = Symbolics.symbolic_linear_solve(expr, p)
@test sol isa Num
@test isequal(sol, y/(y - x))
a, b, islinear = Symbolics.linear_expansion(expr, p)
@test eltype((a, b)) <: Num
@test isequal((a, b, islinear), (-(x - y), -y, true))
@test isequal(unwrap_const(unwrap(Symbolics.symbolic_linear_solve(x * p ~ 0, p))), 0)
@test_throws Any Symbolics.symbolic_linear_solve(1/x + p * p/x ~ 0, p)
@test isequal(Symbolics.symbolic_linear_solve(x * y ~ p, x), p / y)
@test isequal(Symbolics.symbolic_linear_solve(x * -y ~ p, y), -p / x)
@test isequal(Symbolics.symbolic_linear_solve(x * y ~ p, p), x * y)
@test isequal(Symbolics.symbolic_linear_solve(x^2 * y ~ p, y), p / x^2)
@test isequal(Symbolics.symbolic_linear_solve(x^2 * y ~ p, p), x^2 * y)
@test_throws Any Symbolics.symbolic_linear_solve(x^2 * y - sin(p) ~ p, p)
@test Symbolics.symbolic_linear_solve(x^2 * y - sin(p) ~ p, p, check=false) === nothing
@test_throws Any Symbolics.symbolic_linear_solve(t*D(x) ~ y, t)
@test isequal(Symbolics.symbolic_linear_solve(t*D(x) ~ y, D(x)), y/t)
@test isequal(Symbolics.symbolic_linear_solve([t*D(x) ~ y], [D(x)]), [y/t])
@test isequal(Symbolics.symbolic_linear_solve(t*D(x) ~ y, y), t*D(x))
Dx = D(x)
expr = Dx * x + Dx*t - 2//3*x + y*Dx
a, b, islinear = Symbolics.linear_expansion(expr, x)
@test iszero(expand(a * x + b - expr))
@test isequal(Symbolics.symbolic_linear_solve(expr ~ Dx, Dx), (-2//3*x)/(1 - t - x - y))
@test isequal(Symbolics.symbolic_linear_solve(expr ~ Dx, x), (t*Dx + y*Dx - Dx) / ((2//3) - Dx))

exprs = [
 3//2*x + 2y + 10
 7x + 3y - 8
]
xs = [x, y]
A, b, islinear = Symbolics.linear_expansion(exprs, xs)
@test islinear
@test isequal(unwrap_const.(unwrap.(A)), [3//2 2; 7 3])
@test isequal(unwrap_const.(unwrap.(b)), [10; -8])

@variables x y z
eqs = [
        2//1 * x + y - z ~ 2//1
        2//1 + y - z ~ 3//1*x
        2//1 + y - 2z ~ 3//1*z
      ]
@test unwrap_const.(unwrap.([2 1 -1; -3 1 -1; 0 1 -5] * Symbolics.symbolic_linear_solve(eqs, [x, y, z]))) ≈ [2; -2; -2]
@test isequal(Symbolics.symbolic_linear_solve(2//1*x + y - 2//1*z ~ 9//1*x, 1//1*x), (y - 2//1*z) / 7)
@test isequal(Symbolics.symbolic_linear_solve(x + y ~ 0, x), Symbolics.symbolic_linear_solve([x + y ~ 0], x))
@test isequal(Symbolics.symbolic_linear_solve([x + y ~ 0], [x]), Symbolics.symbolic_linear_solve(x + y ~ 0, [x]))
@test isequal(Symbolics.symbolic_linear_solve(2x/z + sin(z), x), sin(z) / (-2 / z))

@testset "linear_expansion with array variables" begin
    @variables x[1:2] y[1:2] z(..)
    @test !Symbolics.linear_expansion(z(x) + x[1], x[1])[3]
    @test !Symbolics.linear_expansion(z(x[1]) + x[1], x[1])[3]
    a, b, islin = Symbolics.linear_expansion(z(x[2]) + x[1], x[1])
    @test islin && isequal(value(a), 1) && isequal(value(b), z(x[2]))
    a, b, islin = Symbolics.linear_expansion((x + x)[1], x[1])
    @test islin && isequal(value(a), 2) && isequal(value(b), 0)
    a, b, islin = Symbolics.linear_expansion(y[1], x[1])
    @test islin && isequal(value(a), 0) && isequal(value(b), y[1])
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

@testset "linear_expansion of matrix-vector products (issue #1231)" begin
    @variables x[1:3] b[1:3] A[1:3, 1:3] s

    Amat, bvec, islin = Symbolics.linear_expansion(A * x + b, x)
    @test islin
    @test isequal(Amat, Symbolics.scalarize(A))
    @test isequal(bvec, Symbolics.scalarize(b))

    Amat, bvec, islin = Symbolics.linear_expansion(A * x + b ~ zeros(Num, 3), x)
    @test islin
    @test isequal(Amat, .-Symbolics.scalarize(A))
    @test isequal(bvec, .-Symbolics.scalarize(b))

    Amat, bvec, islin = Symbolics.linear_expansion(unwrap(A * x + b), unwrap(x))
    @test islin
    @test isequal(Amat, Symbolics.scalarize(A))
    @test isequal(bvec, Symbolics.scalarize(b))

    Amat, bvec, islin = Symbolics.linear_expansion(s .* x .+ b, x)
    @test islin
    @test isequal(Amat, Diagonal(fill(unwrap(s), 3)))
    @test isequal(bvec, Symbolics.scalarize(b))

    # Opaque dependence on the whole array still counts as nonlinear.
    @variables z(..)
    @test !Symbolics.linear_expansion(Symbolics.scalarize(z(x) .+ x), Symbolics.scalarize(x))[3]
end

@testset "linear_expansion array overload preserves LinearExpander semantics" begin
    @variables x y c q[1:3]

    A, b, islin = Symbolics.linear_expansion([x + y], [x + y])
    @test islin
    @test isequal(Num.(A), Num.(reshape([1], 1, 1)))
    @test isequal(Num.(b), Num.([0]))
    @test isequal(Num.(Symbolics.symbolic_linear_solve([0 ~ x + y], [x + y])), Num.([0]))

    A, b, islin = Symbolics.linear_expansion([y * (x + 1)], [x])
    @test islin
    @test isequal(Num.(A), Num.(reshape([y], 1, 1)))
    @test isequal(Num.(b), Num.([y]))
    @test isequal(Num.(Symbolics.symbolic_linear_solve([y * (x + 1) ~ 0], [x])), Num.([-1]))

    A, b, islin = Symbolics.linear_expansion([y * ifelse(c > 0, 2x + 1, 3x + 2)], [x])
    @test islin
    @test isequal(Num.(A), Num.(reshape([y * ifelse(c > 0, 2, 3)], 1, 1)))
    @test isequal(Num.(b), Num.([y * ifelse(c > 0, 1, 2)]))

    A, b, islin = Symbolics.linear_expansion([q[1] * q[2]], [q[1]])
    @test islin
    @test isequal(Num.(A), Num.(reshape([q[2]], 1, 1)))
    @test isequal(Num.(b), Num.([0]))
    @test isequal(Num.(Symbolics.symbolic_linear_solve([q[1] * q[2] ~ 1], [q[1]])), Num.([1 / q[2]]))

    A, b, islin = Symbolics.linear_expansion([y * (3 * (x + y) - 3x)], [x])
    @test islin
    @test isequal(Num.(A), Num.(reshape([0], 1, 1)))
    @test isequal(Num.(b), Num.([3y^2]))

    A, b, islin = Symbolics.linear_expansion([y * (3 * (x + y) - 3x) * (x + 3)], [x])
    @test islin
    @test isequal(Num.(A), Num.(reshape([3y^2], 1, 1)))
    @test isequal(Num.(b), Num.([9y^2]))

    A, b, islin = Symbolics.linear_expansion([(2x + 1) / y], [x])
    @test islin
    @test isequal(Num.(A), Num.(reshape([2 / y], 1, 1)))
    @test isequal(Num.(b), Num.([1 / y]))

    A, b, islin = Symbolics.linear_expansion([sin(y) + 2x], [x])
    @test islin
    @test isequal(Num.(A), Num.(reshape([2], 1, 1)))
    @test isequal(Num.(b), Num.([sin(y)]))

    A, b, islin = Symbolics.linear_expansion([2x + 3y + 4], [x, y])
    @test islin
    @test isequal(Num.(A), Num.([2 3]))
    @test isequal(Num.(b), Num.([4]))

    A, b, islin = Symbolics.linear_expansion([x * y], [x, y])
    @test islin
    @test isequal(Num.(A), Num.([y 0]))
    @test isequal(Num.(b), Num.([0]))
end

@testset "linear_expansion conservative fast path" begin
    @variables x y q[1:2]
    for z0 in (3 * (x + y) - 3x - 3y, 2 * (x + 1) - 2x - 2)
        A, b, islin = Symbolics.linear_expansion([x + x * y * z0, y], [x, y])
        @test islin
        @test isequal(expand.(Num.(A)), Num.([1 0; 0 1]))
        @test isequal(expand.(Num.(b)), Num.([0, 0]))
        sol = Symbolics.symbolic_linear_solve([x + x * y * z0 ~ 1, y ~ 2], [x, y])
        @test isequal(simplify.(sol; expand = true), Num.([1, 2]))
    end
    for (t, xs, expected_A, expected_b) in (
            (x, [x, x], [1 0], [0]),
            (x + y, [x, x + y], [1 0], [y]),
            (y * (x + y), [x, x + y], [y 0], [y^2]),
            (x * y, [x, y], [y 0], [0]),
            (x * y, [y, x], [x 0], [0]),
            (q[1] * q[2], [q[1], q[2]], [q[2] 0], [0]),
        )
        A, b, islin = Symbolics.linear_expansion([t], xs)
        @test islin
        @test isequal(Num.(A), Num.(expected_A))
        @test isequal(Num.(b), Num.(expected_b))
    end
    for t in (x^2, sin(x), y^x)
        @test !Symbolics.linear_expansion([t], [x, y])[3]
    end
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

    @testset "Strict mode" begin
        lex = Symbolics.LinearExpander(p; strict = true)
        @test !lex(ifelse(p < 1, 1, 2))[3]
    end
end

@testset "`linear_expansion` of `LinearAlgebra.dot`" begin
    @variables a b c d
    arr1 = Symbolics.SConst([a, b])
    arr2 = Symbolics.SConst([c, d])
    ex = dot(Symbolics.SConst([arr1, arr2]), [[1, 2], [3, 4]])
    @test isequal(
        Symbolics.scalarize(Symbolics.linear_expansion(ex, a)), 
        (Symbolics.SConst(1), 2b + dot(arr2, [3, 4]), true)
    )
    @test isequal(
        Symbolics.scalarize(Symbolics.linear_expansion(ex, b)), 
        (Symbolics.SConst(2), a + dot(arr2, [3, 4]), true)
    )
    @test isequal(
        Symbolics.scalarize(Symbolics.linear_expansion(ex, c)), 
        (Symbolics.SConst(3), 4d + dot(arr1, [1, 2]), true)
    )
    @test isequal(
        Symbolics.scalarize(Symbolics.linear_expansion(ex, d)),
        (Symbolics.SConst(4), 3c + dot(arr1, [1, 2]), true)
    )
end

@testset "`linear_expansion` of multiplication when multiplicand is singular" begin
    @variables x y

    ex = y * (3 * (x + y) - 3x)
    @test SU.query(isequal(x), unwrap(ex))
    @test isequal(expand(ex), 3y^2)
    a, b, islin = Symbolics.linear_expansion(ex, x)
    @test isequal(a, 0)
    @test isequal(b, 3y^2)
    @test islin

    # Also test case when expression itself is not singular, but a term is
    ex = y * (3 * (x + y) - 3x) * (x + 3)
    @test isequal(expand(ex), 3y^2 * x + 9y^2)
    a, b, islin = Symbolics.linear_expansion(ex, x)
    @test isequal(a, 3y^2)
    @test isequal(b, 9y^2)
    @test islin
end

@testset "`linear_expansion` remainder is free of `x` inside opaque calls" begin
    @variables x h w c g(..)
    ex = g(-h * x + ifelse(c > 0, h * x, h * x) + w)
    a, b, islin = Symbolics.linear_expansion(ex, x)
    @test islin
    @test iszero(a)
    @test isequal(b, g(w))
    @test !SU.query(isequal(unwrap(x)), unwrap(b))
end

matmulwrapper(a, b) = a * b
@register_array_symbolic matmulwrapper(a::AbstractMatrix{Real}, b::AbstractVector{Real}) begin
    size = size(b)
    eltype = eltype(a)
    ndims = 1
end

@testset "`linear_expansion` with `@register_array_symbolic`" begin
    @variables x[1:3] p[1:3, 1:3]
    ex = matmulwrapper(p, x)
    @test !Symbolics.linear_expansion(ex[1], x[1])[3]
end
