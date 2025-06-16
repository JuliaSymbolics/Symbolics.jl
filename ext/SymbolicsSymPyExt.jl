module SymbolicsSymPyExt

if isdefined(Base, :get_extension)
    using Symbolics
    using SymPy
else
    using ..Symbolics
    using ..SymPy
end

using Symbolics: value, symbolics_to_sympy, sympy_to_symbolics
using SymbolicUtils: iscall, operation, arguments, symtype, FnType, Symbolic, Term
using LinearAlgebra

function Symbolics.symbolics_to_sympy(expr)
    expr = value(expr)
    expr isa Symbolic || return expr
    if iscall(expr)
        sop = symbolics_to_sympy(operation(expr))
        sargs = map(symbolics_to_sympy, arguments(expr))
        sop === (^) && length(sargs) == 2 && sargs[2] isa Number ? Base.literal_pow(^, sargs[1], Val(sargs[2])) : sop(sargs...)
    else
        name = string(nameof(expr))
        symtype(expr) <: FnType ? SymPy.SymFunction(name) : SymPy.Sym(name)
    end
end

"""
    sympy_to_symbolics(sympy_expr, vars)
Converts a SymPy expression to Symbolics.jl.
# Arguments
- `sympy_expr`: SymPy expression.
- `vars`: List or dictionary of Symbolics variables.
# Example
```julia
@variables x y
sympy_expr = SymPy.Sym("x")^2 + SymPy.Sym("y")
symbolics_expr = sympy_to_symbolics(sympy_expr, [x, y])
```
"""
function Symbolics.sympy_to_symbolics(sympy_expr, vars)
    if sympy_expr == SymPy.oo
        return Inf
    elseif sympy_expr == -SymPy.oo
        return -Inf
    end
    mod = Module()
    for v in vars
        Core.eval(mod, :($(nameof(v)) = $v))
    end

    Symbolics.parse_expr_to_symbolic(Meta.parse(string(sympy_expr)), mod)
end

"""
    sympy_linear_solve(A, b)
Solves linear system Ax = b using SymPy.
# Arguments
- `A`: Matrix of Symbolics expressions.
- `b`: Vector of Symbolics expressions.
# Returns
Vector of Symbolics solutions.
# Example
```julia
@variables x y
A = [1 2; 3 4]
b = [x, y]
sol = sympy_linear_solve(A, b)
```
"""
function Symbolics.sympy_linear_solve(A, b)
    A_sympy = map(symbolics_to_sympy, A)
    b_sympy = map(symbolics_to_sympy, b)
    A_mat = SymPy.Matrix(A_sympy)
    b_vec = SymPy.Matrix(reshape(b_sympy, length(b_sympy), 1))
    sol_sympy = A_mat \ b_vec
    all_expressions = vcat(vec(A), b)
    vars = Symbolics.get_variables(all_expressions)
    [sympy_to_symbolics(s, vars) for s in sol_sympy]
end

"""
    sympy_algebraic_solve(expr, var)
Solves algebraic equation(s) expr = 0 for var(s) using SymPy.
# Arguments
- `expr`: Symbolics expression or vector of expressions for a system of equations (linear or nonlinear).
- `var`: Symbolics variable or vector of variables to solve for.
# Returns
- For a single equation: List of Symbolics solutions.
- For a system: List of dictionaries mapping variables to solutions.
# Example
```julia
@variables x y
# Single equation
expr = x^2 - 4
sol = sympy_algebraic_solve(expr, x)  # Returns [2, -2]
# Nonlinear system
eqs = [x^2 + y^2 - 4, x - y]  # Circle and line
sol = sympy_algebraic_solve(eqs, [x, y])  # Returns [{x=>1, y=>1}, {x=>-1, y=>-1}]
```
"""
function Symbolics.sympy_algebraic_solve(expr, var)
    expr_sympy = expr isa AbstractVector ? map(symbolics_to_sympy, expr) : symbolics_to_sympy(expr)
    var_sympy = var isa AbstractVector ? map(symbolics_to_sympy, var) : symbolics_to_sympy(var)
    sol_sympy = SymPy.solve(expr_sympy, var_sympy, dict=true)
    vars = var isa AbstractVector ? var : [var]
    if expr isa AbstractVector
        varmap = Dict(string(nameof(v)) => v for v in vars)
        return [Dict(varmap[string(k)] => sympy_to_symbolics(v, vars) for (k, v) in s) for s in sol_sympy]
    else
        return [sympy_to_symbolics(s[var_sympy], vars) for s in sol_sympy]
    end
end

"""
    sympy_integrate(expr, var)
Computes indefinite integral of expr w.r.t. var using SymPy.
# Arguments
- `expr`: Symbolics expression.
- `var`: Symbolics variable.
# Returns
Symbolics integral.
# Example
```julia
@variables x
expr = x^2
result = sympy_integrate(expr, x)
```
"""
function Symbolics.sympy_integrate(expr, var)
    expr_sympy = symbolics_to_sympy(expr)
    var_sympy = symbolics_to_sympy(var)
    result_sympy = SymPy.integrate(expr_sympy, var_sympy)
    sympy_to_symbolics(result_sympy, Symbolics.get_variables(expr))
end

"""
    sympy_limit(expr, var, val)
Computes limit of expr as var approaches val using SymPy.
# Arguments
- `expr`: Symbolics expression.
- `var`: Symbolics variable.
- `val`: Symbolics expression or number.
# Returns
Symbolics limit.
# Example
```julia
@variables x
expr = 1/x
result = sympy_limit(expr, x, 0)
```
"""
function Symbolics.sympy_limit(expr, var, val)
    expr_sympy = symbolics_to_sympy(expr)
    var_sympy = symbolics_to_sympy(var)
    val_sympy = symbolics_to_sympy(val)
    result_sympy = SymPy.limit(expr_sympy, var_sympy => val_sympy)
    sympy_to_symbolics(result_sympy, Symbolics.get_variables(expr))
end

"""
    sympy_simplify(expr)
Simplifies a Symbolics expression using SymPy.
# Arguments
- `expr`: Symbolics expression.
# Returns
Simplified Symbolics expression.
# Example
```julia
@variables x
expr = x^2 + 2x^2
result = sympy_simplify(expr)
```
"""
function Symbolics.sympy_simplify(expr)
    expr_sympy = symbolics_to_sympy(expr)
    result_sympy = SymPy.simplify(expr_sympy)
    sympy_to_symbolics(result_sympy, Symbolics.get_variables(expr))
end

"""
    sympy_ode_solve(expr, func, var)
Solves ODE expr = 0 for function func w.r.t. var using SymPy.
# Arguments
- `expr`: Symbolics expression representing ODE (set to 0).
- `func`: Symbolics function (e.g., f(x)).
- `var`: Independent Symbolics variable.
# Returns
Symbolics solution(s).
# Example
```julia
@variables x
@syms f(x)
expr = Symbolics.Derivative(f, x) - 2*f
sol = sympy_ode_solve(expr, f, x)  # Returns C1*exp(2*x)
```
"""
function Symbolics.sympy_ode_solve(expr, func, var)
    expr_sympy = symbolics_to_sympy(expr)
    func_sympy = symbolics_to_sympy(func)
    var_sympy = symbolics_to_sympy(var)
    sol_sympy = SymPy.dsolve(expr_sympy, func_sympy)
    sympy_to_symbolics(sol_sympy, [func, var])
end

end
