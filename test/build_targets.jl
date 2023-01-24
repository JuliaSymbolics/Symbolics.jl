using Symbolics, Test
using ReferenceTests

@variables t a x(t) y(t)
expr = [a * x - x * y, -3y + x * y]
@test_reference "target_functions/1.stan" Symbolics.build_function(expr, [x, y], [a], t,
                                                                   target = Symbolics.StanTarget())

@test_reference "target_functions/1.c" Symbolics.build_function(expr, [x, y], [a], t,
                                                                target = Symbolics.CTarget(),
                                                                lhsname = :internal_var___du,
                                                                rhsnames = [
                                                                    :internal_var___u,
                                                                    :internal_var___p,
                                                                    :t,
                                                                ])

@test_reference "target_functions/1.m" Symbolics.build_function(expr, [x, y], [a], t,
                                                                target = Symbolics.MATLABTarget())

@test Symbolics.build_function(expr, [x, y], [a], t, target = Symbolics.CTarget(),
                               lhsname = :internal_var___du,
                               rhsnames = [:internal_var___u, :internal_var___p, :t]) ==
      Symbolics.build_function(expr, [x, y], [a], t, target = Symbolics.CTarget(),
                               lhsname = :internal_var___du,
                               rhsnames = [:internal_var___u, :internal_var___p, :t])

@test Symbolics.build_function(expr, [x, y], [a], t, target = Symbolics.StanTarget()) ==
      Symbolics.build_function(expr, [x, y], [a], t, target = Symbolics.StanTarget())

@test Symbolics.build_function(expr, [x, y], [a], t, target = Symbolics.MATLABTarget()) ==
      Symbolics.build_function(expr, [x, y], [a], t, target = Symbolics.MATLABTarget())

@test Symbolics.build_function(expr, [x, y], [a], t, target = Symbolics.CTarget(),
                               lhsname = :internal_var___du,
                               rhsnames = [:internal_var___u, :internal_var___p, :t]) ==
      Symbolics.build_function(expr, vcat(x, y), [a], t, target = Symbolics.CTarget(),
                               lhsname = :internal_var___du,
                               rhsnames = [:internal_var___u, :internal_var___p, :t])

@test Symbolics.build_function(expr, [x, y], [a], t, target = Symbolics.StanTarget()) ==
      Symbolics.build_function(expr, vcat(x, y), [a], t, target = Symbolics.StanTarget())

@test Symbolics.build_function(expr, [x, y], [a], t, target = Symbolics.MATLABTarget()) ==
      Symbolics.build_function(expr, vcat(x, y), [a], t, target = Symbolics.MATLABTarget())

# Matrix CTarget test
let
    @variables x[1:4] y[1:4] z[1:4]
    expression = hcat(x, y, z) # Results in [x1 y1 z1; x2 y2 z2; ...]
    variables = vcat(x, y, z) # Results in [x1, x2, x3, x4, y1, ...]
    cfunc = build_function(expression, variables; target = Symbolics.CTarget(),
                           expression = Val{true})
    @test_reference "target_functions/matrix.c" cfunc
end

# Scalar CTarget test
let
    @variables x y t
    expression = x + y + t
    cfunc = build_function(expression, [x], [y], t; target = Symbolics.CTarget(),
                           expression = Val{true})
    @test_reference "target_functions/scalar1.c" cfunc
end

# Scalar CTarget test with scalar multiplication and powers
let
    @variables x y a t
    expression = x^2 + y^-1 + sin(a)^3.5 + 2t
    cfunc = build_function(expression, [x, y], [a], t; target = Symbolics.CTarget(),
                           expression = Val{true})

    @test_reference "target_functions/scalar2.c" cfunc
end

# Matrix StanTarget test
let
    @variables x[1:4] y[1:4] z[1:4]
    expression = hcat(x, y, z) # Results in [x1 y1 z1; x2 y2 z2; ...]
    variables = vcat(x, y, z) # Results in [x1, x2, x3, x4, y1, ...]
    sfunc = build_function(expression, vcat(x, y), z, []; target = Symbolics.StanTarget(),
                           expression = Val{true})

    @test_reference "target_functions/matrix.stan" sfunc
end

# Scalar StanTarget test
let
    @variables t x(t) y(t) z(t)
    expression = x + y + z
    sfunc = build_function(expression, vcat(x, y), [z], t; target = Symbolics.StanTarget(),
                           expression = Val{true})

    @test_reference "target_functions/scalar1.stan" sfunc
end

# Matrix MATLABTarget test
let
    @variables x[1:4] y[1:4] z[1:4]
    expression = hcat(x, y, z) # Results in [x1 y1 z1; x2 y2 z2; ...]
    variables = vcat(x, y, z) # Results in [x1, x2, x3, x4, y1, ...]
    mfunc = build_function(expression, vcat(x, y, z); target = Symbolics.MATLABTarget(),
                           expression = Val{true})

    @test_reference "target_functions/matrix.m" mfunc
end

# Scalar MATLABTarget test
let
    @variables x y z
    expression = x + y + z
    mfunc = build_function(expression, vcat(x, y, z); target = Symbolics.MATLABTarget(),
                           expression = Val{true})

    @test_reference "target_functions/scalar1.m" mfunc
end
