using Test
using Symbolics
using Symbolics: laplace, inverse_laplace
import Nemo, Groebner

@variables t, s
@syms f(t) F(s)

# https://sites.math.washington.edu/~aloveles/Math307Fall2019/m307LaplacePractice.pdf
@test isequal(laplace(exp(4t) + 5, f, t, s, F), 1/(s-4) + 5/s)
@test isequal(laplace(cos(2t) + 7sin(2t), f, t, s, F), s/(s^2 + 4) + 14/(s^2 + 4))
@test isequal(laplace(exp(-2t)*cos(3t) + 5exp(-2t)*sin(3t), f, t, s, F), (s+2)/((s+2)^2 + 9) + 15/((s+2)^2 + 9))
@test isequal(laplace(10 + 5t + t^2 - 4t^3, f, t, s, F), expand(10/s + 5/s^2 + 2/s^3 - 24/s^4))
@test isequal(laplace(exp(3t)*(t^2 + 4t + 2), f, t, s, F), 2/(s-3)^3 + 4/(s-3)^2 + 2/(s-3))
@test isequal(laplace(6exp(5t)*cos(2t) - exp(7t), f, t, s, F), 6(s-5)/((s-5)^2 + 4) + expand(-1/(s-7)))

# https://www.math.lsu.edu/~adkins/m2065/2065s08review2a.pdf
@test isequal(inverse_laplace(7/(s+3)^3, F, t, s, f), (7//2)t^2 * exp(-3t))
@test_broken isequal(inverse_laplace((s+2)/(s^2 - 3s - 4), F, t, s, f), (7//2)t^2 * exp(-3t)) # partial fraction decomposition