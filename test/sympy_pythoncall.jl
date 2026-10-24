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
@test string(sympy_expr_arr) == "a[1] + a[2]"
back_expr_arr = sympy_pythoncall_to_symbolics(sympy_expr_arr, collect(Symbolics.get_variables(expr_arr)))
@test isequal(Symbolics.simplify(expr_arr), Symbolics.simplify(back_expr_arr))

@variables M[1:2, 1:2]
expr_mat = M[1, 2] + 2M[2, 1]
sympy_expr_mat = symbolics_to_sympy_pythoncall(expr_mat)
@test string(sympy_expr_mat) == "M[1, 2] + 2*M[2, 1]"
back_expr_mat = sympy_pythoncall_to_symbolics(sympy_expr_mat, collect(Symbolics.get_variables(expr_mat)))
@test isequal(Symbolics.simplify(expr_mat), Symbolics.simplify(back_expr_mat))
