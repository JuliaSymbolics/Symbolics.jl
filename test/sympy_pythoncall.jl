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
@test isequal(Symbolics.simplify(result), x^3/3)

# Test 4: Simplification
expr = x^2 + 2x^2
result = sympy_pythoncall_simplify(expr)
@test isequal(Symbolics.simplify(result), 3x^2)

# Test issue #1619: Array element variables map to flat SymPy symbols
@variables a[1:2]
expr_arr = a[1] + a[2]
sympy_expr_arr = symbolics_to_sympy_pythoncall(expr_arr)
@test occursin("a[1]", string(sympy_expr_arr)) && occursin("a[2]", string(sympy_expr_arr))
back_expr_arr = sympy_pythoncall_to_symbolics(sympy_expr_arr, collect(Symbolics.get_variables(expr_arr)))
@test isequal(Symbolics.simplify(expr_arr), Symbolics.simplify(back_expr_arr))
