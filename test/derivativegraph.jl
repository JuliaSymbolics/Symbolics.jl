using Symbolics
using SymbolicUtils
using Test
using Symbolics: value, unwrap, dstar_derivative

@variables x
D = Differential(x)

@test_broken isequal(dstar_derivative(x, x), expand_derivatives(D(x)))
@test isequal(dstar_derivative(2x, x), expand_derivatives(D(2x)))
@test isequal(dstar_derivative(x^2, x), expand_derivatives(D(x^2)))
@test isequal(expand(dstar_derivative(sin(cos(x))*cos(cos(x)), x)), expand_derivatives(D(sin(cos(x))*cos(cos(x)))))
@test isequal(expand(dstar_derivative(cos(sin(exp(2x)) + cos(exp(2x))), x)), expand(expand_derivatives(D(cos(sin(exp(2x)) + cos(exp(2x)))))))
@test isequal(expand(dstar_derivative(2sin(x^4 - x) + 3cos(2x), x)), expand(expand_derivatives(D(2sin(x^4 - x) + 3cos(2x)))))
@test isequal(dstar_derivative(2x*exp(x) + 2*exp(x), x), expand(expand_derivatives(D(2x*exp(x) + 2*exp(x)))))
