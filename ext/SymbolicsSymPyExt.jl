module SymbolicsSymPyExt
if isdefined(Base,:get_extension)
    using Symbolics
    using SymPy
else
    using ..Symbolics
    using ..SymPy
end

using Symbolics: value
using Symbolics.SymbolicUtils
using SymbolicUtils: istree, operation, arguments, symtype,
                        FnType, Symbolic

function Symbolics.symbolics_to_sympy(expr)
    expr = value(expr)
    expr isa Symbolic || return expr
    if istree(expr)
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

end #module
