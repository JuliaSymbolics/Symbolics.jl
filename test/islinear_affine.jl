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

# issue #1105 with a dependent variable: `D(u^2) = 2u*D(u)` is genuinely
# nonlinear in `u`, and querying `t` keeps the conservative handling since
# e.g. `d/dt D(u) = D(D(u))` is nonzero.
let
    @variables t u(t) v(t)
    D = Differential(t)
    @test !Symbolics.isaffine(D(u^2), [u])
    # second derivative of `u * D(u^2)` w.r.t. `u` is nonzero
    @test Matrix(Symbolics.hessian_sparsity(u * D(u^2), [u])) == [true;;]
    @test !Symbolics.isaffine(u * D(u^2), [u])
    @test !Symbolics.islinear(u * D(u^2), [u])
    # second derivative of `t * D(u)` w.r.t. `t` is `2D(D(u)) + t*D(D(D(u)))`
    @test Matrix(Symbolics.hessian_sparsity(t * D(u), [t])) == [true;;]
    @test !Symbolics.isaffine(t * D(u), [t])
    @test !Symbolics.islinear(t * D(u), [t])
    @test Symbolics.isaffine(D(u), [u])
    @test Symbolics.isaffine(D(u) + u, [u])
    @test Symbolics.islinear(u * D(u), [u])
    @test Matrix(Symbolics.hessian_sparsity(D(u * v), [u, v])) == [false true; true false]
end
