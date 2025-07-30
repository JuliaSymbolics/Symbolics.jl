using Test
using SymPyPythonCall
using Symbolics

@variables x y f(x)

# Test 1: Round-trip conversion
expr = x^2 + y
sympy_expr = symbolics_to_sympy_pythoncall(expr)
back_expr = sympy_pythoncall_to_symbolics(sympy_expr, [x, y])
@test isequal(Symbolics.simplify(expr), Symbolics.simplify(back_expr))

# Test 2: Algebraic solver (single equation)
eq = x^2 - 4
sol = sympy_pythoncall_algebraic_solve(eq, x)
@test length(sol) == 2

# Test 3: Integration
expr = x^2
result = sympy_pythoncall_integrate(expr, x)
@test Symbolics.simplify(result - x^3/3) == 0

# Test 4: Simplification
expr = x^2 + 2x^2
result = sympy_pythoncall_simplify(expr)
@test isequal(Symbolics.simplify(result), 3x^2)