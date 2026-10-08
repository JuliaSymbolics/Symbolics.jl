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

# issue #1105: `Differential` applications are opaque to partial
# differentiation (differentiating `D(x)` w.r.t. any variable gives zero),
# so they carry no dependence on the query variables.
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
    @test Symbolics.isaffine(D(x^2), [x])
    @test Symbolics.isaffine(Dx^2, [x])
    @test Symbolics.isaffine(sin(Dx), [x])
end
