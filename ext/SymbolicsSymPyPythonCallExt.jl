module SymbolicsSymPyPythonCallExt

if isdefined(Base, :get_extension)
    using Symbolics
    using SymPyPythonCall
else
    using ..Symbolics
    using ..SymPyPythonCall
end

using Symbolics: value, Differential, Num
using SymbolicUtils: iscall, operation, arguments, symtype, FnType, Symbolic, Term
using LinearAlgebra

# Use the high-level SymPyPythonCall API for conversion
function Symbolics.symbolics_to_sympy_pythoncall(expr)
    expr = value(expr)
    expr isa Symbolic || return expr
    if iscall(expr)
        op = operation(expr)
        args = arguments(expr)

        if op isa Differential
            @assert length(args) == 1 "Differential operator must have exactly one argument."
            return SymPyPythonCall.diff(symbolics_to_sympy_pythoncall(args[1]), symbolics_to_sympy_pythoncall(op.x))
        end

        sop = symbolics_to_sympy_pythoncall(op)
        sargs = map(symbolics_to_sympy_pythoncall, args)
        return sop === (^) && length(sargs) == 2 && sargs[2] isa Number ? sargs[1]^sargs[2] : sop(sargs...)
    else
        name = string(nameof(expr))
        return symtype(expr) <: FnType ? SymPyPythonCall.SymFunction(name) : SymPyPythonCall.Sym(name)
    end
end

function Symbolics.sympy_pythoncall_to_symbolics(sympy_expr, vars)
    # Check for infinity values
    sp = SymPyPythonCall._sympy_
    py_obj = sympy_expr isa SymPyPythonCall.Sym ? sympy_expr.o : sympy_expr
    
    if Bool(SymPyPythonCall.PythonCall.pyeq(py_obj, sp.oo))
        return Inf
    elseif Bool(SymPyPythonCall.PythonCall.pyeq(py_obj, -sp.oo))
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
    
    # Create equations: Ax = b -> Ax - b = 0
    vars_in_system = unique(vcat(Symbolics.get_variables(A), Symbolics.get_variables(b)))
    equations = []
    for i in 1:size(A, 1)
        lhs = sum(A_sympy[i, j] * symbolics_to_sympy_pythoncall(vars_in_system[j]) for j in 1:length(vars_in_system) if j <= size(A, 2))
        eq = lhs - b_sympy[i]
        push!(equations, eq)
    end
    
    # Use high-level SymPyPythonCall solve function
    sol_sympy = SymPyPythonCall.solve(equations, [symbolics_to_sympy_pythoncall(v) for v in vars_in_system])
    all_expressions = vcat(vec(A), b)
    vars = Symbolics.get_variables(all_expressions)
    [sympy_pythoncall_to_symbolics(s, vars) for s in sol_sympy]
end

function Symbolics.sympy_pythoncall_algebraic_solve(expr, var)
    expr_sympy = expr isa AbstractVector ? map(symbolics_to_sympy_pythoncall, expr) : symbolics_to_sympy_pythoncall(expr)
    var_sympy = var isa AbstractVector ? map(symbolics_to_sympy_pythoncall, var) : symbolics_to_sympy_pythoncall(var)
    
    # Use high-level SymPyPythonCall solve function
    sol_sympy = SymPyPythonCall.solve(expr_sympy, var_sympy)
    vars = var isa AbstractVector ? var : [var]
    # Convert each solution
    return [sympy_pythoncall_to_symbolics(s, vars) for s in sol_sympy]
end

function Symbolics.sympy_pythoncall_integrate(expr, var)
    expr_sympy = symbolics_to_sympy_pythoncall(expr)
    var_sympy = symbolics_to_sympy_pythoncall(var)
    # Use high-level SymPyPythonCall integrate function
    result_sympy = SymPyPythonCall.integrate(expr_sympy, var_sympy)
    sympy_pythoncall_to_symbolics(result_sympy, Symbolics.get_variables(expr))
end

function Symbolics.sympy_pythoncall_limit(expr, var, val)
    expr_sympy = symbolics_to_sympy_pythoncall(expr)
    var_sympy = symbolics_to_sympy_pythoncall(var)
    val_sympy = symbolics_to_sympy_pythoncall(val)
    # Use high-level SymPyPythonCall limit function
    result_sympy = SymPyPythonCall.limit(expr_sympy, var_sympy, val_sympy)
    sympy_pythoncall_to_symbolics(result_sympy, Symbolics.get_variables(expr))
end

function Symbolics.sympy_pythoncall_simplify(expr)
    expr_sympy = symbolics_to_sympy_pythoncall(expr)
    # Use high-level SymPyPythonCall simplify function
    result_sympy = SymPyPythonCall.simplify(expr_sympy)
    sympy_pythoncall_to_symbolics(result_sympy, Symbolics.get_variables(expr))
end

function Symbolics.sympy_pythoncall_ode_solve(expr, func, var)
    expr_sympy = symbolics_to_sympy_pythoncall(expr)
    func_sympy = symbolics_to_sympy_pythoncall(func)
    var_sympy = symbolics_to_sympy_pythoncall(var)
    # Use high-level SymPyPythonCall dsolve function
    sol_sympy = SymPyPythonCall.dsolve(expr_sympy, func_sympy)
    sol_expr = SymPyPythonCall.rhs(sol_sympy)
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