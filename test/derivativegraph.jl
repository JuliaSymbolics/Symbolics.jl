using Symbolics
using SymbolicUtils
using Test
using Symbolics: value, unwrap

@variables x y

test_roots = [x, 2x, x^2, sin(cos(x))*cos(cos(x)), cos(sin(exp(2x)) + cos(exp(2x))), 2sin(x^4 - x) + 3cos(2x), 2x*exp(x) + 2*exp(x)]

dgs = []

for root in test_roots
    dg = Symbolics.DerivativeGraph([unwrap(root)], [unwrap(x)], Int)
    push!(dgs, dg)
end