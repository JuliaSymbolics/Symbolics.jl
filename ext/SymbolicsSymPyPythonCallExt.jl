module SymbolicsSymPyPythonCallExt

if isdefined(Base, :get_extension)
    using Symbolics
    using SymPyPythonCall
else
    using ..Symbolics
    using ..SymPyPythonCall
end

using Symbolics: value, symbolics_to_sympy_pythoncall, sympy_pythoncall_to_symbolics, Differential, Num
using SymbolicUtils: iscall, operation, arguments, symtype, FnType, Symbolic, Term
using LinearAlgebra

# Access the underlying SymPy Python module
const sp = SymPyPythonCall._sympy_
const PythonCall = SymPyPythonCall.PythonCall

function Symbolics.symbolics_to_sympy_pythoncall(expr)
    expr = value(expr)
    expr isa Symbolic || return expr
    if iscall(expr)
        op = operation(expr)
        args = arguments(expr)

        if op isa Differential
            @assert length(args) == 1 "Differential operator must have exactly one argument."
            return sp.Derivative(symbolics_to_sympy_pythoncall(args[1]), symbolics_to_sympy_pythoncall(op.x))
        end

        sop = symbolics_to_sympy_pythoncall(op)
        sargs = map(symbolics_to_sympy_pythoncall, args)
        return sop === (^) && length(sargs) == 2 && sargs[2] isa Number ? Base.literal_pow(^, sargs[1], Val(sargs[2])) : sop(sargs...)
    else
        name = string(nameof(expr))
        return symtype(expr) <: FnType ? sp.Function(name) : sp.Symbol(name)
    end
end

function Symbolics.sympy_pythoncall_to_symbolics(sympy_expr, vars)
    # Check for infinity values using pyeq and convert to Julia boolean
    if Bool(PythonCall.pyeq(sympy_expr, sp.oo))
        return Inf
    elseif Bool(PythonCall.pyeq(sympy_expr, -sp.oo))
        return -Inf
    end
    dict = Dict{Symbol, Any}()
    for v in vars
        dict[nameof(v)] = v
    end

    # Convert Python/SymPy notation to Julia notation
    expr_str = string(sympy_expr)
    expr_str = replace(expr_str, "**" => "^")  # Convert exponentiation
    Symbolics.parse_expr_to_symbolic(Meta.parse(expr_str), dict)
end

function Symbolics.sympy_pythoncall_linear_solve(A, b)
    A_sympy = map(symbolics_to_sympy_pythoncall, A)
    b_sympy = map(symbolics_to_sympy_pythoncall, b)
    
    # Convert A to Python list of lists
    A_py_rows = [PythonCall.pylist(row) for row in eachrow(A_sympy)]
    A_py = PythonCall.pylist(A_py_rows)
    A_mat = sp.Matrix(A_py)
    
    # Convert b to Python list and then column matrix
    b_py = PythonCall.pylist(b_sympy)
    b_vec = sp.Matrix(length(b_sympy), 1, b_py)
    
    sol_sympy = A_mat.LUsolve(b_vec)
    all_expressions = vcat(vec(A), b)
    vars = Symbolics.get_variables(all_expressions)
    [sympy_pythoncall_to_symbolics(s, vars) for s in sol_sympy]
end

function Symbolics.sympy_pythoncall_algebraic_solve(expr, var)
    expr_sympy = expr isa AbstractVector ? map(symbolics_to_sympy_pythoncall, expr) : symbolics_to_sympy_pythoncall(expr)
    var_sympy = var isa AbstractVector ? map(symbolics_to_sympy_pythoncall, var) : symbolics_to_sympy_pythoncall(var)
    
    # Convert arrays to Python lists for SymPy
    if expr isa AbstractVector
        expr_sympy = PythonCall.pylist(expr_sympy)
    end
    if var isa AbstractVector
        var_sympy = PythonCall.pylist(var_sympy)
    end
    
    sol_sympy = sp.solve(expr_sympy, var_sympy, dict=true)
    vars = var isa AbstractVector ? var : [var]
    if expr isa AbstractVector
        varmap = Dict(string(nameof(v)) => v for v in vars)
        return [Dict(varmap[string(k)] => sympy_pythoncall_to_symbolics(v, vars) for (k, v) in s) for s in sol_sympy]
    else
        return [sympy_pythoncall_to_symbolics(s[var_sympy], vars) for s in sol_sympy]
    end
end

function Symbolics.sympy_pythoncall_integrate(expr, var)
    expr_sympy = symbolics_to_sympy_pythoncall(expr)
    var_sympy = symbolics_to_sympy_pythoncall(var)
    result_sympy = sp.integrate(expr_sympy, var_sympy)
    sympy_pythoncall_to_symbolics(result_sympy, Symbolics.get_variables(expr))
end

function Symbolics.sympy_pythoncall_limit(expr, var, val)
    expr_sympy = symbolics_to_sympy_pythoncall(expr)
    var_sympy = symbolics_to_sympy_pythoncall(var)
    val_sympy = symbolics_to_sympy_pythoncall(val)
    result_sympy = sp.limit(expr_sympy, var_sympy, val_sympy)
    sympy_pythoncall_to_symbolics(result_sympy, Symbolics.get_variables(expr))
end

function Symbolics.sympy_pythoncall_simplify(expr)
    expr_sympy = symbolics_to_sympy_pythoncall(expr)
    result_sympy = sp.simplify(expr_sympy)
    sympy_pythoncall_to_symbolics(result_sympy, Symbolics.get_variables(expr))
end

function Symbolics.sympy_pythoncall_ode_solve(expr, func, var)
    expr_sympy = symbolics_to_sympy_pythoncall(expr)
    func_sympy = symbolics_to_sympy_pythoncall(func)
    var_sympy = symbolics_to_sympy_pythoncall(var)
    sol_sympy = sp.dsolve(expr_sympy, func_sympy)
    sol_expr = sol_sympy.rhs
    parsing_vars = Vector{SymbolicUtils.BasicSymbolic}()
    vars_in_expr = Symbolics.get_variables(value(expr))
    func_val = value(func)
    for v in vars_in_expr
        if !isequal(v, func_val)
            push!(parsing_vars, v)
        end
    end
    push!(parsing_vars, value(var))
    push!(parsing_vars, operation(func_val))
    unwrapped_vars = unique(parsing_vars)
    sympy_pythoncall_to_symbolics(sol_expr, unwrapped_vars)
end

end