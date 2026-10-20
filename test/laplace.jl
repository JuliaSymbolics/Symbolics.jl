using Test
using Symbolics
using Symbolics: laplace, inverse_laplace
import Nemo, Groebner

@variables t, s
@syms f(t) F(s)

# https://sites.math.washington.edu/~aloveles/Math307Fall2019/m307LaplacePractice.pdf
@test isequal(laplace(exp(4t) + 5, f, t, F, s), 1/(s-4) + 5/s)
@test isequal(laplace(cos(2t) + 7sin(2t), f, t, F, s), s/(s^2 + 4) + 14/(s^2 + 4))
@test isequal(laplace(exp(-2t)*cos(3t) + 5exp(-2t)*sin(3t), f, t, F, s), (s+2)/((s+2)^2 + 9) + 15/((s+2)^2 + 9))
@test isequal(laplace(10 + 5t + t^2 - 4t^3, f, t, F, s), expand(10/s + 5/s^2 + 2/s^3 - 24/s^4))
@test isequal(laplace(exp(3t)*(t^2 + 4t + 2), f, t, F, s), 2/(s-3)^3 + 4/(s-3)^2 + 2/(s-3))
@test isequal(laplace(6exp(5t)*cos(2t) - exp(7t), f, t, F, s), 6(s-5)/((s-5)^2 + 4) + expand(-1/(s-7)))

# regression tests for fixed transform rules
@test isequal(laplace(t*exp(t), f, t, F, s), 1/(s-1)^2)
@test isequal(laplace(t^3*exp(t), f, t, F, s), 6/(s-1)^4)
@test isequal(laplace(t*sin(t), f, t, F, s), 2s/(1 + s^2)^2)
@test isequal(laplace(sqrt(t), f, t, F, s), Symbolics.term(sqrt, pi)/(2 * s^(3//2)))
@test isequal(laplace(t^(1//2), f, t, F, s), Symbolics.term(sqrt, pi)/(2 * s^(3//2)))
@test isequal(laplace(t^(3//2), f, t, F, s), Symbolics.term(Symbolics.gamma, 5//2)/s^(5//2))
@test isequal(laplace(t^(5//2), f, t, F, s), Symbolics.term(Symbolics.gamma, 7//2)/s^(7//2))
@test isequal(laplace(t^(-1//2), f, t, F, s), Symbolics.term(Symbolics.gamma, 1//2)/s^(1//2))
# powers with n <= -1 diverge: fall back to the integral definition rather than erroring
@test Symbolics._has_integral(laplace(t^(-2), f, t, F, s))
@test Symbolics._has_integral(laplace(1/t, f, t, F, s))
@test isequal(laplace(t^(1//2)*exp(2t), f, t, F, s), Symbolics.term(sqrt, pi)/(2 * (s - 2)^(3//2)))
@test isequal(laplace(t^(3//2)*exp(2t), f, t, F, s), Symbolics.term(Symbolics.gamma, 5//2)/(s - 2)^(5//2))
# frequency shift with a coefficient-free exp(t)
@test isequal(laplace(t^(1//2)*exp(t), f, t, F, s), Symbolics.term(sqrt, pi)/(2 * (s - 1)^(3//2)))
@test isequal(laplace(t^(3//2)*exp(t), f, t, F, s), Symbolics.term(Symbolics.gamma, 5//2)/(s - 1)^(5//2))
@test isequal(laplace(exp(t)*sin(t), f, t, F, s), 1/((s - 1)^2 + 1))
@test isequal(laplace(exp(t)*cos(2t), f, t, F, s), (s - 1)/((s - 1)^2 + 4))
@test isequal(laplace(exp(t)*f(t), f, t, F, s), F(s - 1))
# negative powers are stored as quotients, e.g. exp(2t)/t^(1//2)
@test isequal(laplace(t^(-1//2)*exp(2t), f, t, F, s), Symbolics.term(Symbolics.gamma, 1//2)/(s - 2)^(1//2))
@test isequal(laplace(t^(-1//2)*exp(t), f, t, F, s), Symbolics.term(Symbolics.gamma, 1//2)/(s - 1)^(1//2))
@test Symbolics._has_integral(laplace(t^(-3//2)*exp(2t), f, t, F, s))
# s-derivative rule
@test isequal(laplace(t^2*f(t), f, t, F, s), (Differential(s)^2)(F(s)))
# fractional or symbolic powers without a closed-form rule fall back to the integral definition
@test Symbolics._has_integral(laplace(t^(1//2)*f(t), f, t, F, s))
@variables n
@test Symbolics._has_integral(laplace(t^n, f, t, F, s))

# opt-in trig_reduce for products and powers of trig/hyperbolic functions (values verified against SymPy)
@test isequal(laplace(2cos(t)*sin(t), f, t, F, s; trig_reduce=true), 2/(s^2 + 4))
@test isequal(laplace(sin(t)^2, f, t, F, s; trig_reduce=true), -(1//2)*s/(s^2 + 4) + (1//2)/s)
@test isequal(laplace(cos(t)^3, f, t, F, s; trig_reduce=true), (1//4)*s/(s^2 + 9) + (3//4)*s/(s^2 + 1))
@test isequal(laplace(sin(t)*cos(2t), f, t, F, s; trig_reduce=true), -(1//2)/(s^2 + 1) + (3//2)/(s^2 + 9))
@test isequal(laplace(sin(2t)*sin(3t), f, t, F, s; trig_reduce=true), -(1//2)*s/(s^2 + 25) + (1//2)*s/(s^2 + 1))
@test isequal(laplace(sinh(t)*cosh(t), f, t, F, s; trig_reduce=true), 1/(s^2 - 4))
@test isequal(laplace(3 + sin(t)^2, f, t, F, s; trig_reduce=true), -(1//2)*s/(s^2 + 4) + (1//2)/s + 3/s)
@test isequal(laplace(exp(2t)*sin(t)^2, f, t, F, s; trig_reduce=true), -(1//2)*(s - 2)/((s - 2)^2 + 4) + (1//2)/(s - 2))
# off by default; when on, a term trig_reduce cannot make transformable keeps the original fallback
@test Symbolics._has_integral(laplace(2cos(t)*sin(t), f, t, F, s))
@test isequal(laplace(tan(t)^2, f, t, F, s; trig_reduce=true), laplace(tan(t)^2, f, t, F, s))

# https://www.math.lsu.edu/~adkins/m2065/2065s08review2a.pdf
@test isequal(inverse_laplace(7/(s+3)^3, F, s, f, t), (7//2)t^2 * exp(-3t))
@test isequal(inverse_laplace((s-9)/(s^2 + 9), F, s, f, t), cos(3t) - 3sin(3t))
# regression test for fixed inverse rule
@test isequal(inverse_laplace(s*(-4 + s^2)/(4 + s^2)^2, F, s, f, t; warns=false), cos(2t) - 2t*sin(2t))
# partial fraction decomposition
@test isequal(inverse_laplace((s+2)/(s^2 - 3s - 4), F, s, f, t), (6//5)*exp(4t) - (1//5)*exp(-t))
@test isequal(inverse_laplace(1/(s^2 - 10s + 9), F, s, f, t), (1//8)*exp(9t) - (1//8)*exp(t))

Dt = Differential(t)
@test_skip isequal(laplace_solve_ode(Dt(f(t)) + 3f(t) ~ t^2*exp(-3t) + t*exp(-2t) + t, f, t, [1]), (1//3)*t^3*exp(-3t) + t*exp(-2t) + (1//3)*t + (19//9)*exp(-3t) - exp(-2t) - 1//9) # generates very large system of equations that solve_interms_ofvar has trouble with
@test isequal(laplace_solve_ode((Dt^2)(f(t)) - 3Dt(f(t)) + 2f(t) ~ 4, f, t, [2, 3]), 2 - 3exp(t) + 3exp(2t))
@test isequal(laplace_solve_ode((Dt^3)(f(t)) - Dt(f(t)) ~ 2, f, t, [4,4,4]), 5exp(t) - exp(-t) - 2t)
@test isequal(laplace_solve_ode((Dt^3)(f(t)) - Dt(f(t)) ~ 6 - 3t^2, f, t, [1, 1, 1]), exp(t) + t^3)
@test isequal(laplace_solve_ode((Dt^2)(f(t)) - f(t) ~ 2sin(t), f, t, [0, 0]), (1//2)exp(t) - (1//2)exp(-t) - sin(t))
@test isequal(laplace_solve_ode((Dt^2)(f(t)) + 2Dt(f(t)) ~ 5f(t), f, t, [0, 0]), 0)
@test isequal(laplace_solve_ode((Dt^2)(f(t)) + f(t) ~ sin(4t), f, t, [0, 0]), (4//15)sin(t) - (1//15)sin(4t))
@test isequal(laplace_solve_ode((Dt^2)(f(t)) + Dt(f(t)) ~ 1 + 2t, f, t, [0, 0]), 1 - exp(-t) + t^2 - t)
@test isequal(laplace_solve_ode((Dt^2)(f(t)) + 4Dt(f(t)) + 3f(t) ~ 6, f, t, [0, 0]), exp(-3t) - 3exp(-t) + 2)
@test isequal(laplace_solve_ode((Dt^2)(f(t)) - 2Dt(f(t)) ~ 3*(t + exp(2t)), f, t, [0, 0]), (3//8) - (3//4)t - (3//4)t^2 - (3//8)exp(2t) + (3//2)t*exp(2t))
@test_broken isequal(laplace_solve_ode((Dt^2)(f(t)) - 2Dt(f(t)) ~ 20*exp(-t)*cos(t), f, t, [0, 0], warns=false), 3exp(2t) - 5 + 2exp(-t)*cos(t) - 4exp(-t)*sin(t)) # irreducible quadratic in inverse laplace
@test isequal(laplace_solve_ode((Dt^2)(f(t)) + f(t) ~ 2 + 2cos(t), f, t, [0, 0]), 2 - 2cos(t) + t*sin(t))
@test isequal(laplace_solve_ode((Dt^2)(f(t)) - Dt(f(t)) ~ 30cos(3t), f, t, [0, 0]), 3exp(t) - 3cos(3t) - sin(3t))