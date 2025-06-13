module SymbolicsSymPyExt

if isdefined(Base, :get_extension)
    using Symbolics
    using SymPy
else
    using ..Symbolics
    using ..SymPy
end

using Symbolics: value
using SymbolicUtils: iscall, operation, arguments, symtype, FnType, Symbolic
using LinearAlgebra

# Existing symbolics_to_sympy function
function Symbolics.symbolics_to_sympy(expr)
    expr = value(expr)
    expr isa Symbolic || return expr
    if iscall(expr)
        sop = symbolics_to_sympy(operation(expr))
        sargs = map(symbolics_to_sympy, arguments(expr))
        if sop === (^) && length(sargs) == 2 && sargs[2] isa Number
            return Base.literal_pow(^, sargs[1], Val(sargs[2]))
        else
            return sop(sargs...)
        end
    else # isa Symbolics.Sym
        name = string(nameof(expr))
        return symtype(expr) <: FnType ? SymPy.SymFunction(name) : SymPy.Sym(name)
    end
end

"""
    sympy_to_symbolics(sympy_expr, vars)

Converts a SymPy expression to a Symbolics.jl expression.

This function takes a SymPy expression and a list or dictionary of Symbolics variables
to ensure proper variable mapping during conversion. It is used to return solver results
to the Symbolics.jl format.

# Arguments
- `sympy_expr`: The SymPy expression to convert (typically a `SymPy.Sym` object).
- `vars`: A list or dictionary mapping variable names to Symbolics.jl variables.

# Example
```julia
@variables x y
sympy_expr = SymPy.Sym("x")^2 + SymPy.Sym("y")
symbolics_expr = sympy_to_symbolics(sympy_expr, [x, y])
```

See also: [`symbolics_to_sympy`]
"""
function sympy_to_symbolics(sympy_expr, vars)
    varmap = Dict(string(nameof(v)) => v for v in vars)
    # Convert SymPy expression to string and parse back to Symbolics
    expr_str = string(sympy_expr)
    # Use Symbolics' parsing with variable mapping
    return Symbolics.parse(expr_str, varmap)
end

# Wrapper Functions for SymPy Solvers

"""
    sympy_linear_solve(A, b)

Solves a linear system Ax = b using SymPy's linear solver.

# Arguments
- `A`: A matrix of Symbolics.jl expressions.
- `b`: A vector of Symbolics.jl expressions.

# Returns
A vector of Symbolics.jl expressions representing the solution.

# Example
```julia
@variables x y
A = [1 2; 3 4]
b = [x, y]
sol = sympy_linear_solve(A, b)
```
"""
function sympy_linear_solve(A, b)
    A_sympy = map(symbolics_to_sympy, A)
    b_sympy = map(symbolics_to_sympy, b)
    A_mat = SymPy.SymMatrix(A_sympy)
    b_vec = SymPy.SymMatrix(b_sympy)
    sol_sympy = A_mat \ b_vec  # SymPy's linear solver
    vars = unique(vcat(vec(Symbolics.get_variables.(A)), Symbolics.get_variables.(b)))
    return [sympy_to_symbolics(s, vars) for s in sol_sympy]
end

"""
    sympy_algebraic_solve(expr, var)

Solves an algebraic equation expr = 0 for the variable var using SymPy.

# Arguments
- `expr`: A Symbolics.jl expression representing the equation.
- `var`: The Symbolics.jl variable to solve for.

# Returns
A list of Symbolics.jl expressions representing the solutions.

# Example
```julia
@variables x
expr = x^2 - 4
sol = sympy_algebraic_solve(expr, x)
```
"""
function sympy_algebraic_solve(expr, var)
    expr_sympy = symbolics_to_sympy(expr)
    var_sympy = symbolics_to_sympy(var)
    sol_sympy = SymPy.solve(expr_sympy, var_sympy)
    return [sympy_to_symbolics(s, [var]) for s in sol_sympy]
end

"""
    sympy_integrate(expr, var)

Computes the indefinite integral of expr with respect to var using SymPy.

# Arguments
- `expr`: A Symbolics.jl expression to integrate.
- `var`: The Symbolics.jl variable of integration.

# Returns
A Symbolics.jl expression representing the integral.

# Example
```julia
@variables x
expr = x^2
result = sympy_integrate(expr, x)
```
"""
function sympy_integrate(expr, var)
    expr_sympy = symbolics_to_sympy(expr)
    var_sympy = symbolics_to_sympy(var)
    result_sympy = SymPy.integrate(expr_sympy, var_sympy)
    vars = Symbolics.get_variables(expr)
    return sympy_to_symbolics(result_sympy, vars)
end

"""
    sympy_limit(expr, var, val)

Computes the limit of expr as var approaches val using SymPy.

# Arguments
- `expr`: A Symbolics.jl expression.
- `var`: The Symbolics.jl variable.
- `val`: The value the variable approaches (Symbolics or numeric).

# Returns
A Symbolics.jl expression or value representing the limit.

# Example
```julia
@variables x
expr = 1/x
result = sympy_limit(expr, x, 0)
```
"""
function sympy_limit(expr, var, val)
    expr_sympy = symbolics_to_sympy(expr)
    var_sympy = symbolics_to_sympy(var)
    val_sympy = symbolics_to_sympy(val)
    result_sympy = SymPy.limit(expr_sympy, var_sympy, val_sympy)
    vars = Symbolics.get_variables(expr)
    return sympy_to_symbolics(result_sympy, vars)
end

"""
    sympy_simplify(expr)

Simplifies a Symbolics.jl expression using SymPy's simplification.

# Arguments
- `expr`: A Symbolics.jl expression to simplify.

# Returns
A Symbolics.jl expression representing the simplified result.

# Example
```julia
@variables x
expr = x^2 + 2x^2
result = sympy_simplify(expr)
```
"""
function sympy_simplify(expr)
    expr_sympy = symbolics_to_sympy(expr)
    result_sympy = SymPy.simplify(expr_sympy)
    vars = Symbolics.get_variables(expr)
    return sympy_to_symbolics(result_sympy, vars)
end

# Tests
"""
    run_sympy_tests()

Runs tests to verify Symbolics <-> SymPy round-trip conversion and solver wrappers.
"""
function run_sympy_tests()
    @variables x y

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
    @test length(sol) == 2
    @test all(s -> isequal(s, 2) || isequal(s, -2), sol)

    # Test 4: Integration
    expr = x^2
    result = sympy_integrate(expr, x)
    @test isequal(result, x^3/3)

    # Test 5: Limit
    expr = 1/x
    result = sympy_limit(expr, x, 0)
    @test isequal(result, Symbolics.Infinity())

    # Test 6: Simplification
    expr = x^2 + 2x^2
    result = sympy_simplify(expr)
    @test isequal(result, 3x^2)

    println("All SymPy extension tests passed!")
end

# Exports
export symbolics_to_sympy, sympy_to_symbolics
export sympy_linear_solve, sympy_algebraic_solve, sympy_integrate, sympy_limit, sympy_simplify
export run_sympy_tests

end # module