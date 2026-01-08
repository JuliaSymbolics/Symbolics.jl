module SymbolicsLatexifyExt

using Symbolics
using Latexify
using LaTeXStrings
using TermInterface
using SymbolicUtils
using Symbolics: value, hide_lhs, postwalk, wrap
using SymbolicUtils: BSImpl, FnType, unwrap, symtype, BasicSymbolic
using Moshi.Match: @match

# metadata to specify how to format syms
struct SymLatexWrapper end
Symbolics.option_to_metadata_type(::Val{:latexwrapper}) = SymLatexWrapper

function default_latex_wrapper(sym)
    length(sym) <= 1 && return sym
    return string("\\mathtt{", sym, "}")
end

prettify_expr(expr) = expr
prettify_expr(f::Function) = nameof(f)
prettify_expr(expr::Expr) = Expr(expr.head, prettify_expr.(expr.args)...)

function cleanup_exprs(ex)
    return postwalk(x -> iscall(x) && length(arguments(x)) == 0 ? operation(x) : x, ex)
end

function latexify_derivatives(ex)
    return postwalk(ex) do x
        Meta.isexpr(x, :call) || return x
        if x.args[1] == :_derivative
            num, den, deg = x.args[2:end]
            if num isa Expr && length(num.args) == 2
                return Expr(:call, :/,
                            Expr(:call, :*,
                                 "\\mathrm{d}$(deg == 1 ? "" : "^{$deg}")", num
                                ),
                            diffdenom(den)
                           )
            else
                return Expr(:call, :*,
                            Expr(:call, :/,
                                 "\\mathrm{d}$(deg == 1 ? "" : "^{$deg}")",
                                 diffdenom(den)
                                ),
                            num
                           )
            end
    elseif x.args[1] === :_integral
        lower, upper, var_of_int, integrand = x.args[2:end]
        lower_s = strip(latexify(lower).s, '\$')
        upper_s = strip(latexify(upper).s, '\$')
        return Expr(:call, :*,
                "\\int_{$lower_s}^{$upper_s)}",
                var_of_int,
                integrand
            )
        elseif x.args[1] === :_textbf
            ls = latexify(latexify_derivatives(sorted_arguments(x)[1])).s
            return "\\textbf{" * strip(ls, '\$') * "}"
        else
            return x
        end
    end
end

recipe(n) = latexify_derivatives(cleanup_exprs(_toexpr(n)))

@latexrecipe function f(n::Num)
    env --> :equation
    mult_symbol --> "~"
    fmt --> FancyNumberFormatter(5)
    index --> :subscript
    snakecase --> true
    safescripts --> true

    return recipe(value(n))
end

@latexrecipe function f(z::Complex{Num})
    env --> :equation
    mult_symbol --> "~"
    index --> :subscript

    iszero(z.im) && return :($(recipe(value(z.re))))
    iszero(z.re) && return :($(recipe(value(z.im))) * $im)
    return :($(recipe(value(z.re))) + $(recipe(value(z.im))) * $im)
end

@latexrecipe function f(n::Function)
    env --> :equation
    mult_symbol --> "~"
    index --> :subscript

    return nameof(n)
end


@latexrecipe function f(n::Symbolics.Arr)
    env --> :equation
    mult_symbol --> "~"
    index --> :subscript

    return value(n)
end

@latexrecipe function f(n::Symbolics.CallAndWrap)
    env --> :equation
    mult_symbol --> "~"
    index --> :subscript

    return n.f
end

@latexrecipe function f(n::SymbolicUtils.BasicSymbolic)
    env --> :equation
    mult_symbol --> "~"
    index --> :subscript

    return recipe(n)
end

@latexrecipe function f(eqs::Vector{Equation})
    index --> :subscript
    has_connections = any(x -> hide_lhs(value(x.lhs)), eqs)
    if has_connections
        env --> :equation
        return map(first∘first∘Latexify.apply_recipe, eqs)
    else
        env --> :align
        return wrap.(getfield.(eqs, :lhs)), wrap.(getfield.(eqs, :rhs))
    end
end

@latexrecipe function f(eq::Equation)
    env --> :equation
    index --> :subscript

    if hide_lhs(value(eq.lhs)) || !(value(eq.lhs) isa Union{Number, AbstractArray, BasicSymbolic})
        return value(eq.rhs)
    else
        return Expr(:(=), recipe(eq.lhs), recipe(eq.rhs))
    end
end

Base.show(io::IO, ::MIME"text/latex", x::Symbolics.RCNum) = print(io, "\$\$ " * latexify(x) * " \$\$")
Base.show(io::IO, ::MIME"text/latex", x::SymbolicUtils.BasicSymbolic) = print(io, "\$\$ " * latexify(x) * " \$\$")
Base.show(io::IO, ::MIME"text/latex", x::Equation) = print(io, "\$\$ " * latexify(x) * " \$\$")
Base.show(io::IO, ::MIME"text/latex", x::Vector{Equation}) = print(io, "\$\$ " * latexify(x) * " \$\$")
Base.show(io::IO, ::MIME"text/latex", x::AbstractArray{<:Symbolics.RCNum}) = print(io, "\$\$ " * latexify(x) * " \$\$")

# `_toexpr` is only used for latexify
function _toexpr(O; latexwrapper = default_latex_wrapper)
    O = unwrap(O)
    SymbolicUtils.isconst(O) && return value(O)
    if SymbolicUtils.ismul(O)
        m = O
        numer = Any[]
        denom = Any[]

        # We need to iterate over each term in m, ignoring the numeric coefficient.
        # This iteration needs to be stable, so we can't iterate over m.dict.
        for term in Iterators.drop(sorted_arguments(m), isone(m.coeff) ? 0 : 1)
            if !SymbolicUtils.ispow(term)
                push!(numer, _toexpr(term))
                continue
            end
            base, pow = arguments(term)
            pow = value(pow)
            isneg = (pow isa Number && pow < 0) || (iscall(pow) && operation(pow) === (-) && length(arguments(pow)) == 1)
            if !isneg
                if SymbolicUtils._isone(pow)
                    pushfirst!(numer, _toexpr(base))
                else
                    pushfirst!(numer, Expr(:call, :^, _toexpr(base), _toexpr(pow)))
                end
            else
                newpow = -1*pow
                if SymbolicUtils._isone(newpow)
                    pushfirst!(denom, _toexpr(base))
                else
                    pushfirst!(denom, Expr(:call, :^, _toexpr(base), _toexpr(newpow)))
                end
            end
        end

        if !isreal(m.coeff)
            numer_expr = Expr(:call, :*, m.coeff, numer...)
        elseif isempty(numer) || !isone(abs(m.coeff))
            numer_expr = Expr(:call, :*, abs(m.coeff), numer...)
        else
            numer_expr = length(numer) > 1 ? Expr(:call, :*, numer...) : numer[1]
        end

        if isempty(denom)
            frac_expr = numer_expr
        else
            denom_expr = length(denom) > 1 ? Expr(:call, :*, denom...) : denom[1]
            frac_expr = Expr(:call, :/, numer_expr, denom_expr)
        end

        if isreal(m.coeff) && real(m.coeff) < 0
            return Expr(:call, :-, frac_expr)
        else
            return frac_expr
        end
    end
    if SymbolicUtils.issym(O) 
        sym = string(nameof(O))
        sym = replace(sym, Symbolics.NAMESPACE_SEPARATOR => ".")

        # override if the sym has its own latex wrapper
        symwrapper = hasmetadata(O, SymLatexWrapper) ? getmetadata(O, SymLatexWrapper) : 
            latexwrapper
        sym = symwrapper(sym)
        return Symbol(sym)
    end
    !iscall(O) && return O

    op = operation(O)
    args = sorted_arguments(O)
    latexwrapper = hasmetadata(O, SymLatexWrapper) ? getmetadata(O, SymLatexWrapper) : 
        default_latex_wrapper

    if (op===(*)) && (args[1] === -1)
        arg_mul = Expr(:call, :(*), _toexpr(args[2:end])...)
        return Expr(:call, :(-), arg_mul)
    end

    if op isa Differential
      #  DERIVATIVES LOGIC
        num = args[1]
        diff_var = op.x 
        
        deg = op.order
        
        while iscall(num) && operation(num) isa Differential && isequal(operation(num).x, diff_var)
            inner_op = operation(num)
            deg += inner_op.order
            num = arguments(num)[1]
        end

        den = deg > 1 ? (diff_var^deg) : diff_var
        return :(_derivative($(_toexpr(num)), $den, $deg))   

    elseif op isa Integral
        lower = op.domain.domain.left
        upper = op.domain.domain.right
        vars = op.domain.variables
        integrand = args[1]
        var = if vars isa Tuple
            Expr(:call, :(*), _toexpr(vars...))
        else
            _toexpr(vars)
        end
        return Expr(:call, :_integral, _toexpr(lower), _toexpr(upper), vars, _toexpr(integrand))
    elseif symtype(op) <: FnType
        isempty(args) && return nameof(op)
        return Expr(:call, _toexpr(op; latexwrapper), _toexpr(args)...)
    elseif op === getindex && symtype(args[1]) <: AbstractArray
        return getindex_to_symbol(O)
    elseif op === (\)
        return :(solve($(_toexpr(args[1])), $(_toexpr(args[2]))))
    elseif SymbolicUtils.issym(op) && SymbolicUtils.symtype(op) <: AbstractArray
        return :(_textbf($(nameof(op))))
    elseif op === identity
        return _toexpr(only(args)) # suppress identity transformations (e.g. "identity(π)" -> "π")
    end
    return Expr(:call, Symbol(op), _toexpr(args; latexwrapper)...)
end
_toexpr(x::Integer; latexwrapper = default_latex_wrapper) = x
_toexpr(x::AbstractFloat; latexwrapper = default_latex_wrapper) = x

function _toexpr(eq::Equation; latexwrapper = default_latex_wrapper)
    Expr(:(=), _toexpr(eq.lhs), _toexpr(eq.rhs))
end

_toexpr(eqs::AbstractArray; latexwrapper = default_latex_wrapper) = map(eq->_toexpr(eq), eqs)
_toexpr(x::Num; latexwrapper = default_latex_wrapper) = _toexpr(value(x))

function getindex_to_symbol(t)
    @assert iscall(t) && operation(t) === getindex && SymbolicUtils.symtype(sorted_arguments(t)[1]) <: AbstractArray
    args = sorted_arguments(t)
    idxs = args[2:end]
    O = args[1]
    latexwrapper = hasmetadata(O, SymLatexWrapper) ? getmetadata(O, SymLatexWrapper) : 
        default_latex_wrapper

    # this is to ensure X(t)[1] becomes X_1(t) in Latex
    if iscall(O) && SymbolicUtils.issym(operation(O))
        oop = operation(O)        
        oargs = sorted_arguments(O)
        return :($(_toexpr(oop; latexwrapper))[$(idxs...)]($(_toexpr(oargs)...)))
    else
        return :($(_toexpr(O; latexwrapper))[$(idxs...)])
    end
end

function diffdenom(e)
    e = unwrap(e)
    if SymbolicUtils.issym(e)
        LaTeXString("\\mathrm{d}$e")
    elseif SymbolicUtils.ispow(e)
        base, expo = arguments(e)
        suffix = SymbolicUtils._isone(expo) ? "" : "^{$(expo)}"
        LaTeXString("\\mathrm{d}$(base)$(suffix)")
   elseif SymbolicUtils.ismul(e)
        LaTeXString(prod(diffdenom(arg).s for arg in arguments(e)))
    else
        LaTeXString("\\mathrm{d}$e")
    end
end

end
