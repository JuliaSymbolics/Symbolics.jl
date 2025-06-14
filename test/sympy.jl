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
using SymbolicsSymPyExt

@variables t p(t) x y(t)
expr = x * p + (x^2 - 1 + y) * (p + 2t)
sexpr = symbolics_to_sympy(expr)
sp = symbolics_to_sympy(p)

@test SymPy.simplify(symbolics_to_sympy(Symbolics.solve_for(expr, p))) == SymPy.solve(sexpr, sp)[1]

@variables x y
@syms f(x)

# Test 1: Round-trip conversion
expr = x^2 + y
sympy_expr = symbolics_to_sympy(expr)
back_expr = sympy_to_symbolics(sympy_expr, [x, y])
@test isequal(expr, back_expr)

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
eqs = [x^2 + y^2 - 4, x - y]  # Circle and line
sol_sys = sympy_algebraic_solve(eqs, [x, y])
@test length(sol_sys) == 2 && any(d -> isequal(d[x], 1) && isequal(d[y], 1), sol_sys)

# Test 5: Integration
expr = x^2
result = sympy_integrate(expr, x)
@test isequal(result, x^3/3)

# Test 6: Limit
expr = 1/x
result = sympy_limit(expr, x, 0)
@test isequal(result, Symbolics.Infinity())

# Test 7: Simplification
expr = x^2 + 2x^2
result = sympy_simplify(expr)
@test isequal(result, 3x^2)

# Test 8: ODE solver
D = Differential(x)
ode = D(f) - 2*f
sol_ode = sympy_ode_solve(ode, f, x)
@test isequal(sol_ode, Symbolics.parse("C1*exp(2*x)", Dict("f"=>f, "x"=>x)))