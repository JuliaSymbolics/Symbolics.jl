using Symbolics
using Symbolics: solve_linear_ode, LinearODE, has_const_coeffs, to_homogeneous, symbolic_solve_ode, find_particular_solution, IVP, solve_IVP
using Groebner, Nemo
using Test

@variables x, y, t
Dt = Symbolics.Differential(t)

# Systems
# ideally, `isapprox` would all be `isequal`, but there seem to be some floating point inaccuracies
@test isapprox(solve_linear_ode([1 0; 0 -1], [1, -1], t), [exp(t), -exp(-t)])
@test isapprox(solve_linear_ode([-3 4; -2 3], [7, 2], t), [10exp(-t) - 3exp(t), 5exp(-t) - 3exp(t)])
@test isapprox(solve_linear_ode([4 -3; 8 -6], [7, 2], t), [18 - 11exp(-2t), 24 - 22exp(-2t)])

@test_broken isapprox(solve_linear_ode([-1 -2; 2 -1], [1, -1], t), [exp(-t)*(cos(2t) + sin(2t)), exp(-t)*(sin(2t) - cos(2t))]) # can't handle complex eigenvalues (though it should be able to)

@test isapprox(solve_linear_ode([1 -1 0; 1 2 1; -2 1 -1], [7, 2, 3], t), (5//3)*exp(-t)*[-1, -2, 7] - 14exp(t)*[-1, 0, 1] + (16//3)*exp(2t)*[-1, 1, 1])

@test isequal(solve_linear_ode([1 0; 0 -1], [1, -1], t), [exp(t), -exp(-t)])
@test isequal(solve_linear_ode([-3 4; -2 3], [7, 2], t), [10exp(-t) - 3exp(t), 5exp(-t) - 3exp(t)])
@test isapprox(solve_linear_ode([4 -3; 8 -6], [7, 2], t), [18 - 11exp(-2t), 24 - 22exp(-2t)])

@test isequal(solve_linear_ode([1 -1 0; 1 2 1; -2 1 -1], [7, 2, 3], t), (5//3)*exp(-t)*[-1, -2, 7] - 14exp(t)*[-1, 0, 1] + (16//3)*exp(2t)*[-1, 1, 1])

@test_throws ArgumentError solve_linear_ode([1 2; 3 4], [1, 2, 3], t) # mismatch between A and x0
@test_throws ArgumentError solve_linear_ode([1 2 3; 4 5 6], [1, 2], t) # A isn't square
@test_throws ArgumentError Symbolics.solve_uncoupled_system([1 2; 3 4], [1, 2], t) # A isn't diagonal

# LinearODEs
@test Symbolics.is_homogeneous(LinearODE(x, t, [1, 1], 0))
@test !Symbolics.is_homogeneous(LinearODE(x, t, [t, 1], t^2))

@test has_const_coeffs(LinearODE(x, t, [1, 1], 0))
@test !has_const_coeffs(LinearODE(x, t, [t^2, 1], 0))

@test Symbolics.is_homogeneous(to_homogeneous(LinearODE(x, t, [t, 1], t^2)))
@test !Symbolics.is_linear_ode(((Dt^2)(x))^2 ~ x^3, x, t)

C = Symbolics.variables(:C, 1:5)

## constant coefficients, nth-order
@test isequal(symbolic_solve_ode(LinearODE(x, t, [-1], 0)), C[1]*exp(t))
@test isequal(symbolic_solve_ode(LinearODE(x, t, [-4, 3], 0)), C[1]*exp(-4t) + C[2]*exp(t))

## first order (solving via integrating factor can be found in test/sympy.jl)
@test isequal(symbolic_solve_ode(LinearODE(x, t, [1], 2sin(t))), C[1]*exp(-t) + sin(t) - cos(t))

## repeated characteristic roots
@test isequal(symbolic_solve_ode((Dt^2)(x) + 2(Dt^1)(x) + x ~ 0, x, t), C[1]*exp(-t) + C[2]*t*exp(-t))
@test isequal(symbolic_solve_ode(LinearODE(x, t, [0, 0, 0, 4, -4], 0)), C[1] + C[2]*t + C[3]*t^2 + C[4]*exp(2t) + C[5]*t*exp(2t))
@test isequal(symbolic_solve_ode(LinearODE(x, t, [8, 12, 6], 0)), C[1]*exp(-2t) + C[2]*t*exp(-2t) + C[3]*t^2*exp(-2t))

## complex characteristic roots
@test isequal(symbolic_solve_ode(LinearODE(x, t, [1, 0], 0)), C[1]*cos(t) + C[2]*sin(t))
@test isequal(symbolic_solve_ode(LinearODE(x, t, [0, 1, 0], 0)), C[1] + C[2]*cos(t) + C[3]*sin(t))

## resonant response formula
@test isequal(symbolic_solve_ode(LinearODE(x, t, [9, -6], 4exp(3t))), C[1]*exp(3t) + C[2]*t*exp(3t) + 2(t^2)*exp(3t))
### trig functions
@test isequal(symbolic_solve_ode(LinearODE(x, t, [6, 5], 2exp(-t)*cos(t))), C[1]*exp(-2t) + C[2]*exp(-3t) + (1//5)*exp(-t)*cos(t)+(3//5)*exp(-t)*sin(t))

## undetermined coefficients
@test isequal(symbolic_solve_ode(LinearODE(x, t, [-3, 2], 2t - 5)), C[1]exp(t) + C[2]exp(-3t) - (2//3)t + 11//9)
@test isequal(find_particular_solution(LinearODE(x, t, [1, 0], t^2)), t^2 - 2)

# Parsing
@test isequal(LinearODE(x, t, [1], 0), LinearODE(Dt(x) + x ~ 0, x, t))
@test isequal(LinearODE(x, t, [sin(t), 0, 3t^2], exp(2t) + 2cos(t)), LinearODE(6t^2*(Dt^2)(x) + 2sin(t)*x - 2exp(2t) + 2(Dt^3)(x) ~ 4cos(t), x, t))

# IVP
@test isequal(solve_IVP(IVP(LinearODE(x, t, [-3, 2], 0), [1, -1])), (1//2)exp(-3t) + (1//2)exp(t))
@test isequal(solve_IVP(IVP(LinearODE(x, t, [9, -6], 4exp(3t)), [5, 6])), 5exp(3t) - 9t*exp(3t) + 2(t^2)*exp(3t))

# Other methods

## Clairaut's equation
@test isequal(symbolic_solve_ode(x ~ Dt(x)*t - ((Dt(x))^3), x, t), C[1]*t - C[1]^3)
@test isequal(symbolic_solve_ode(Dt(x)*t + (Dt(x))^2 - sin(Dt(x)) + 2 ~ x, x, t), C[1]*t + C[1]^2 - sin(C[1]) + 2)
@test isnothing(symbolic_solve_ode(Dt(x) + (Dt(x))^2 ~ x, x, t))
@test isnothing(symbolic_solve_ode(Dt(x)*t + 2t*(Dt(x))^2 ~ x, x, t))
@test isnothing(symbolic_solve_ode(Dt(x) + x*(Dt(x))^2 ~ x, x, t))

## Bernoulli equations (integrating factor solve in test/sympy.jl)
@test isequal(symbolic_solve_ode(Dt(x) - 5x ~ exp(-2t)*x^(-2), x, t), (C[1]exp(15t) - (3//17)exp(-2t))^(1//3))
@test isnothing(symbolic_solve_ode(sqrt((Dt^4)(x)) ~ log(x)^t, x, t))

# Helper function tests
ys = Symbolics.variables(:y, 1:2)
@test isequal(Symbolics.reduce_order((Dt^2)(x) + 3Dt(x) + 2x ~ 0, x, t, ys), [ys[2], -2ys[1] - 3ys[2]])
@test isequal(Symbolics.unreduce_order([ys[1], ys[2]], x, t, ys), [x, Dt(x)])

@test Symbolics.is_solution(C[1]*exp(3t) + C[2]*t*exp(3t) + 2(t^2)*exp(3t), LinearODE(x, t, [9, -6], 4exp(3t)))
@test Symbolics.is_solution(C[1]*exp(-t) + C[2]*t*exp(-t), (Dt^2)(x) + 2(Dt^1)(x) + x ~ 0, x, t)