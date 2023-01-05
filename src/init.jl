function __init__()
    @require SymPy="24249f21-da20-56a4-8eb1-6a02cf4ae2e6" begin

        using Symbolics
        using Symbolics: value
        using SymbolicUtils: istree, operation, arguments, symtype,
                             FnType, Symbolic
        function symbolics_to_sympy(expr)
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

    end # SymPy

    specialize_methods()
end
