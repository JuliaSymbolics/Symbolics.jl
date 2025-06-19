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

function Symbolics.sympy_to_symbolics(sympy_expr, vars)
    if sympy_expr == SymPy.oo
        return Inf
    elseif sympy_expr == -SymPy.oo
        return -Inf
    end
    dict = Dict{Symbol, Any}()
    for v in vars
        dict[nameof(v)] = v
    end

    Symbolics.parse_expr_to_symbolic(Meta.parse(string(sympy_expr)), dict)
end

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

function Symbolics.sympy_integrate(expr, var)
    expr_sympy = symbolics_to_sympy(expr)
    var_sympy = symbolics_to_sympy(var)
    result_sympy = SymPy.integrate(expr_sympy, var_sympy)
    sympy_to_symbolics(result_sympy, Symbolics.get_variables(expr))
end

function Symbolics.sympy_limit(expr, var, val)
    expr_sympy = symbolics_to_sympy(expr)
    var_sympy = symbolics_to_sympy(var)
    val_sympy = symbolics_to_sympy(val)
    result_sympy = SymPy.limit(expr_sympy, var_sympy => val_sympy)
    sympy_to_symbolics(result_sympy, Symbolics.get_variables(expr))
end

function Symbolics.sympy_simplify(expr)
    expr_sympy = symbolics_to_sympy(expr)
    result_sympy = SymPy.simplify(expr_sympy)
    sympy_to_symbolics(result_sympy, Symbolics.get_variables(expr))
end

function Symbolics.sympy_ode_solve(expr, func, var)
    expr_sympy = symbolics_to_sympy(expr)
    func_sympy = symbolics_to_sympy(func)
    var_sympy = symbolics_to_sympy(var)
    sol_sympy = SymPy.dsolve(expr_sympy, func_sympy)
    sympy_to_symbolics(sol_sympy, [func, var])
end

end
