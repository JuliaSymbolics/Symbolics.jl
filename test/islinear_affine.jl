using Symbolics, Test
@variables t, x(t), y(t), z(t)
@test Symbolics.islinear(x + y,[x,y])
@test Symbolics.islinear(x,[x,y])
@test Symbolics.islinear(y,[x,y])
@test !Symbolics.islinear(z,[x,y])

@test Symbolics.isaffine(x + y,[x,y])
@test Symbolics.isaffine(x,[x,y])
@test Symbolics.isaffine(y,[x,y])
@test Symbolics.isaffine(z,[x,y])
@test Symbolics.isaffine(x + y + z,[x,y])

@test Symbolics.isaffine(x + z * y,[x,y])
@test Symbolics.islinear(x + z * y,[x,y])
@test Symbolics.islinear(z * x + z * y,[x,y])
@test Symbolics.islinear(z * (x + y),[x,y])
@test Symbolics.isaffine(z * (x + y),[x,y])

@test Symbolics.isaffine(ifelse(x < 1, y, z), [z])
@test !Symbolics.isaffine(ifelse(x < 1, x, z), [x])

@test !Symbolics.isaffine(ifelse(y < 2, ifelse(x < 1, 0, 1), y), [x])

# issue #1105
let
    @variables x y t
    D = Differential(t)
    Dx = D(x)
    expr = Dx * x + Dx * t - 2 // 3 * x + y * Dx
    @test Symbolics.isaffine(expr - Dx, [x])
    # nonzero constant term `Dx * (t + y - 1)`
    @test !Symbolics.islinear(expr - Dx, [x])
    @test isempty(Symbolics.hessian_sparsity(expr - Dx, [x]).nzval)
    @test Symbolics.isaffine(Dx, [x])
    @test !Symbolics.islinear(Dx, [x])
    @test Symbolics.isaffine(Dx + x, [x])
    @test Symbolics.isaffine(x * Dx, [x])
    @test Symbolics.islinear(x * Dx, [x])
    # nonlinear argument: `D` keeps the argument's dependence
    @test !Symbolics.isaffine(D(x^2), [x])
    @test Symbolics.isaffine(Dx^2, [x])
    @test Symbolics.isaffine(sin(Dx), [x])
end

# expected patterns: Hessian of `expand_derivatives(ex)` with `D(u)`, `D(v)`, ... as independent coordinates
let
    @variables t x u(t) v(t) a(t)
    D = Differential(t)
    covers(S, H) = size(S) == size(H) && all(Matrix(S) .>= H)
    cases = [
        (x * D(t * x), [x], [true;;]),
        (D(t * x)^2, [x], [true;;]),
        (u * D(t * u), [u], [true;;]),
        (u * D(v * u), [u], [true;;]),
        (u * D(u * v), [u], [true;;]),
        (u * D(a * u), [u], [true;;]),
        (D(v * u)^2, [u], [true;;]),
        (u * D(D(u) * u), [u], [true;;]),
        (sin(D(t * u)), [u], [true;;]),
        (v * D(v * u), [u, v], [false true; true true]),
        (u * D(u^2), [u], [true;;]),
        (t * D(u), [t], [true;;]),
    ]
    for (ex, vars, H) in cases
        @test covers(Symbolics.hessian_sparsity(ex, vars), H)
        @test !Symbolics.isaffine(ex, vars)
        @test !Symbolics.islinear(ex, vars)
    end
    @test Symbolics.isaffine(x * D(2x + 3) + D(D(x))^2, [x])
    @test Symbolics.isaffine(x * D(x / 2 - x + t), [x])
end
