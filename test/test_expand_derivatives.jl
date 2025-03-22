using Test
using Symbolics

# Define variables
@variables a x
D = Differential(x)

# Test 1: Rational instead of Float64
ex = exp(a * x) / (2a)
result = expand_derivatives(D(D(ex)))
@test result == (1//2) * a * exp(a * x)

# Test 2: Ensure existing expressions remain unaffected
simple_ex = a * x^2
result_simple = expand_derivatives(D(D(simple_ex)))
@test result_simple == 2a

# Test 3: Ensure float-to-rational conversion works
ex_float = exp(a * x) / 2.0
result_float = expand_derivatives(D(D(ex_float)))
@test result_float == (1//2) * a * exp(a * x)

println("All tests passed!")
