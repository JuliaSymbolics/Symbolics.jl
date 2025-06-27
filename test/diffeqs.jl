using Symbolics
using Symbolics: firstorder_separable_ode_solve, solve_linear_system, LinearODE, is_homogeneous, has_const_coeffs, to_homogeneous, homogeneous_solve
import Groebner, Nemo, SymPy
using Test

@variables x, y, t, C
Dt = Differential(t)
@test_broken isequal(firstorder_separable_ode_solve(Dt(x) - x ~ 0, x, t), C*exp(t))
@test isequal(firstorder_separable_ode_solve(Dt(x) + 2*t*x ~ 0, x, t), C*exp(-(t^2)))
@test isequal(firstorder_separable_ode_solve(Dt(x) + (t^2-3)*x ~ 0, x, t), C*exp(3t - (1//3)*t^3))

# Systems
@test isapprox(solve_linear_system([1 0; 0 -1], [1, -1], t), [exp(t), -exp(-t)])
@test isapprox(solve_linear_system([-3 4; -2 3], [7, 2], t), [10exp(-t) - 3exp(t), 5exp(-t) - 3exp(t)])
@test isapprox(solve_linear_system([4 -3; 8 -6], [7, 2], t), [18 - 11exp(-2t), 24 - 22exp(-2t)])

@test_broken isapprox(solve_linear_system([-1 -2; 2 -1], [1, -1], t), [exp(-t)*(cos(2t) + sin(2t)), exp(-t)*(sin(2t) - cos(2t))])

@test isapprox(solve_linear_system([1 -1 0; 1 2 1; -2 1 -1], [7, 2, 3], t), (5//3)*exp(-t)*[-1, -2, 7] - 14exp(t)*[-1, 0, 1] + (16//3)*exp(2t)*[-1, 1, 1])

@test isequal(solve_linear_system([1 0; 0 -1], [1, -1], t), [exp(t), -exp(-t)])
@test isequal(solve_linear_system([-3 4; -2 3], [7, 2], t), [10exp(-t) - 3exp(t), 5exp(-t) - 3exp(t)])
@test_broken isequal(solve_linear_system([4 -3; 8 -6], [7, 2], t), [18 - 11exp(-2t), 24 - 22exp(-2t)])

@test_broken isequal(solve_linear_system([-1 -2; 2 -1], [1, -1], t), [exp(-t)*(cos(2t) + sin(2t)), exp(-t)*(sin(2t) - cos(2t))])

@test isequal(solve_linear_system([1 -1 0; 1 2 1; -2 1 -1], [7, 2, 3], t), (5//3)*exp(-t)*[-1, -2, 7] - 14exp(t)*[-1, 0, 1] + (16//3)*exp(2t)*[-1, 1, 1])

# LinearODEs
@test is_homogeneous(LinearODE(x, t, [1, 1], 0))
@test !is_homogeneous(LinearODE(x, t, [t, 1], t^2))

@test has_const_coeffs(LinearODE(x, t, [1, 1], 0))
@test !has_const_coeffs(LinearODE(x, t, [t^2, 1], 0))

@test is_homogeneous(to_homogeneous(LinearODE(x, t, [t, 1], t^2)))
@variables C[1:3]
@test isequal(homogeneous_solve(LinearODE(x, t, [-1], 0)), C[1]*exp(t))
@test isequal(homogeneous_solve(LinearODE(x, t, [-4, 3], 0)), C[1]*exp(-4t) + C[2]*exp(t))