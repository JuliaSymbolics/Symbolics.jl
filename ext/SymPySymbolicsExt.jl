module SymPySymbolics
if isdefined(Base,:get_extension)
    using Symbolics
    using PythonCall
    using CondaPkg
else
    using ..Symbolics
    using ..PythonCall
    using ..CondaPkg
end

# rule functions
function pyconvert_rule_sympy_symbol(::Type{Symbolics.Num}, x::Py)
    if !pyisinstance(x,sp.Symbol)
        return PythonCall.pyconvert_unconverted()
    end
    name = PythonCall.pyconvert(Symbol,x.name)
    return PythonCall.pyconvert_return(Symbolics.variable(name))
end

function pyconvert_rule_sympy_pow(::Type{Symbolics.Num}, x::Py)
    if !pyisinstance(x,sp.Pow)
        return PythonCall.pyconvert_unconverted()
    end
    expbase = pyconvert(Symbolics.Num,x.base)
    exp = pyconvert(Symbolics.Num,x.exp)
    return PythonCall.pyconvert_return(expbase^exp)
end

function pyconvert_rule_sympy_mul(::Type{Symbolics.Num}, x::Py)
    if !pyisinstance(x,sp.Mul)
        return PythonCall.pyconvert_unconverted()
    end
    mult = reduce(*,PythonCall.pyconvert.(Symbolics.Num,x.args))
    return PythonCall.pyconvert_return(mult)
end

function pyconvert_rule_sympy_add(::Type{Symbolics.Num}, x::Py)
    if !pyisinstance(x,sp.Add)
       return PythonCall.pyconvert_unconverted() 
    end
    sum = reduce(+, PythonCall.pyconvert.(Symbolics.Num,x.args))
    return PythonCall.pyconvert_return(sum)
end

function pyconvert_rule_sympy_equality(::Type{Symbolics.Equation}, x::Py)
    if !pyisinstance(x,sp.Equality)
        return PythonCall.pyconvert_unconverted()
    end
    rhs = pyconvert(Symbolics.Num,x.rhs)
    lhs = pyconvert(Symbolics.Num,x.lhs)
    return PythonCall.pyconvert_return(rhs ~ lhs)
end

function pyconvert_rule_sympy_derivative(::Type{Symbolics.Num}, x::Py)
    if !pyisinstance(x,sp.Derivative)
        return PythonCall.pyconvert_unconverted()
    end
    variables = pyconvert.(Symbolics.Num,x.variables)
    derivatives = prod(var -> Differential(var), variables)
    expr = pyconvert(Symbolics.Num, x.expr)
    return PythonCall.pyconvert_return(derivatives(expr))
end

function pyconvert_rule_sympy_function(::Type{Symbolics.Num}, x::Py)
    if !pyisinstance(x,sp.Function)
        return PythonCall.pyconvert_unconverted()
    end
    name = pyconvert(Symbol,x.name)
    args = pyconvert.(Symbolics.Num,x.args)
    func = @variables $name(..)
    return PythonCall.pyconvert_return(first(func)(args...))
end

# added rules
PythonCall.pyconvert_add_rule("sympy.core.power:Pow", Symbolics.Num, pyconvert_rule_sympy_pow)

PythonCall.pyconvert_add_rule("sympy.core.symbol:Symbol", Symbolics.Num, pyconvert_rule_sympy_symbol)

PythonCall.pyconvert_add_rule("sympy.core.mul:Mul", Symbolics.Num, pyconvert_rule_sympy_mul)

PythonCall.pyconvert_add_rule("sympy.core.add:Add", Symbolics.Num, pyconvert_rule_sympy_add)

PythonCall.pyconvert_add_rule("sympy.core.relational:Equality", Symbolics.Equation, pyconvert_rule_sympy_equality)

PythonCall.pyconvert_add_rule("sympy.core.function:Derivative",Symbolics.Num, pyconvert_rule_sympy_derivative)

PythonCall.pyconvert_add_rule("sympy.core.function:Function",Symbolics.Num, pyconvert_rule_sympy_function)

function sympy_to_symbolics(expr::Py)
    if pyconvert(Bool,pytype(expr) == sp.core.relational.Equality)
        return pyconvert(Symbolics.Equation,expr)
    end
    return pyconvert(Symbolics.Num,expr)
end

end