if get(ENV, "CI", nothing) !== nothing
    # copied and modified from https://github.com/tkf/IPython.jl/blob/master/test/install_dependencies.jl

    # Adding Pkg in test/REQUIRE would be an error in 0.6.  Using
    # Project.toml still has some gotchas.  So:
    Pkg = Base.require(Base.PkgId(Base.UUID(0x44cfe95a1eb252eab672e2afdf69b78f), "Pkg"))

    # Let PyCall.jl use Python interpreter from Conda.jl
    # See: https://github.com/JuliaPy/PyCall.jl
    ENV["PYTHON"] = ""
    Pkg.build("PyCall")
end

using Test

using SymPy
using Symbolics

@variables t p(t) x y(t)
expr = x * p + (x^2 - 1 + y) * (p + 2t)
sexpr = symbolics_to_sympy(expr)
sp = symbolics_to_sympy(p)

@test SymPy.simplify(symbolics_to_sympy(Symbolics.symbolic_linear_solve(expr, p))) ==
      SymPy.simplify(SymPy.solve(sexpr, sp)[1])

symbolics_sol = SymPy.simplify(symbolics_to_sympy(Symbolics.symbolic_linear_solve(expr, p)))
sympy_sols = SymPy.solve(SymPy.expand(sexpr), sp)
@test !isempty(sympy_sols) && isequal(symbolics_sol, sympy_sols[1])

@variables x y f(x)

# Test 1: Round-trip conversion
expr = x^2 + y
sympy_expr = symbolics_to_sympy(expr)
back_expr = sympy_to_symbolics(sympy_expr, [x, y])
@test isequal(Symbolics.simplify(expr), Symbolics.simplify(back_expr))

# Test 2: Linear solver
A = [1 2; 3 4]
b = [x, y]
sol = sympy_linear_solve(A, b)
@test length(sol) == 2

# Test 3: Algebraic solver
eq = x^2 - 4
sol = sympy_algebraic_solve(eq, x)
@test length(sol) == 2 && all(s -> isequal(s, 2) || isequal(s, -2), sol)

# Test 4: System of equations
eqs = [x^2 + y^2 - 4, x - y]
sol_sys = sympy_algebraic_solve(eqs, [x, y])
@test length(sol_sys) == 2 && any(
    d -> isapprox(Symbolics.substitute(d[x], Dict()), sqrt(2)) &&
         isapprox(Symbolics.substitute(d[y], Dict()), sqrt(2)),
    sol_sys)

# Test 5: Integration
expr = x^2
result = sympy_integrate(expr, x)
@test isequal(Symbolics.simplify(result), x^3/3)

# Test 6: Limit
expr = 1/x
result = sympy_limit(expr, x, 0)
@test isequal(result, Inf)

# Test 7: Simplification
expr = x^2 + 2x^2
result = sympy_simplify(expr)
@test isequal(Symbolics.simplify(result), 3x^2)

# Test 8: ODE solver
D = Symbolics.Differential(x)
@variables C1
ode = D(f) - 2*f
sol_ode = sympy_ode_solve(ode, f, x)
sol_vars = Symbolics.get_variables(sol_ode)
const_sym = only(filter(v -> startswith(string(Symbolics.nameof(v)), "C"), sol_vars))
expected_sol = C1 * exp(2 * x)
canonical_sol_ode = Symbolics.substitute(sol_ode, Dict(const_sym => C1))
@test isequal(canonical_sol_ode, expected_sol)

## Native ODE solver, but using sympy_integrate
@variables x, t
Dt = Symbolics.Differential(t)
C = Symbolics.variables(:C, 1:5)

@test isequal(symbolic_solve_ode(SymbolicLinearODE(x, t, [5/t], 7t)), Symbolics.sympy_simplify(C[1]*t^(-5) + t^2))
@test isequal(symbolic_solve_ode(SymbolicLinearODE(x, t, [cos(t)], cos(t))), 1 + C[1]*exp(-sin(t)))
@test isequal(symbolic_solve_ode(SymbolicLinearODE(x, t, [-(1+t)], 1+t)), Symbolics.expand(Symbolics.sympy_simplify(C[1]*exp((1//2)t^2 + t) - 1)))
# SymPy is being weird and not simplifying correctly (and some symbols are wrong, like pi and erf being syms), but these otherwise work
@test_broken isequal(symbolic_solve_ode(SymbolicLinearODE(x, t, [-2t], 1)), Symbolics.sympy_simplify(exp(t^2)*sqrt(Symbolics.variable(:pi))*erf(t)/2 + C[1]*exp(t^2)))

## Bernoulli equations
@test isequal(symbolic_solve_ode(Dt(x) + (4//t)*x ~ t^3 * x^2, x, t), 1/(C[1]t^4 - t^4 * log(t)))
@test isnothing(symbolic_solve_ode(Dt(x) + sqrt(t)*x ~ sin(4t)*x^3, x, t))

# Test issue #1605: Power with symbolic exponent
@variables x y
expr_power = x^y
sympy_expr_power = symbolics_to_sympy(expr_power)
back_expr_power = sympy_to_symbolics(sympy_expr_power, [x, y])
@test isequal(Symbolics.simplify(expr_power), Symbolics.simplify(back_expr_power))
