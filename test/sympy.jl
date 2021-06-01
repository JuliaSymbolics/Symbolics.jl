using Test
using SymPy
using Symbolics

@variables t p(t) x y(t)
expr = x * p + (x^2 - 1 + y) * (p + 2t)
sexpr = symbolics_to_sympy(expr)
sp = symbolics_to_sympy(p)

@test SymPy.simplify(symbolics_to_sympy(Symbolics.solve_for(expr, p))) == SymPy.solve(sexpr, sp)[1]
