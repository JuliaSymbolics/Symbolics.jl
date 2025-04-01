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
    cdot --> false
    fmt --> FancyNumberFormatter(5)
    index --> :subscript
    snakecase --> true
    safescripts --> true

    return recipe(n)
end

@latexrecipe function f(z::Complex{Num})
    env --> :equation
    cdot --> false
    index --> :subscript

    iszero(z.im) && return :($(recipe(z.re)))
    iszero(z.re) && return :($(recipe(z.im)) * $im)
    return :($(recipe(z.re)) + $(recipe(z.im)) * $im)
end

@latexrecipe function f(n::ArrayOp)
    env --> :equation
    cdot --> false
    index --> :subscript
    return recipe(n.term)
end

@latexrecipe function f(n::Function)
    env --> :equation
    cdot --> false
    index --> :subscript

    return nameof(n)
end


@latexrecipe function f(n::Arr)
    env --> :equation
    cdot --> false
    index --> :subscript

    return unwrap(n)
end

@latexrecipe function f(n::Symbolic)
    env --> :equation
    cdot --> false
    index --> :subscript

    return recipe(n)
end

@latexrecipe function f(eqs::Vector{Equation})
    index --> :subscript
    has_connections = any(x->x.lhs isa Connection, eqs)
    if has_connections
        env --> :equation
        return map(first∘first∘Latexify.apply_recipe, eqs)
    else
        env --> :align
        return Num.(getfield.(eqs, :lhs)), Num.(getfield.(eqs, :rhs))
    end
end

@latexrecipe function f(eq::Equation)
    env --> :equation
    index --> :subscript

    if !(eq.lhs isa Union{Number, AbstractArray, Symbolic})
        return eq.rhs
    else
        return Expr(:(=), Num(eq.lhs), Num(eq.rhs))
    end
end

@latexrecipe function f(c::Connection)
    index --> :subscript
    return Expr(:call, :connect, map(nameof, c.systems)...)
end

Base.show(io::IO, ::MIME"text/latex", x::RCNum) = print(io, "\$\$ " * latexify(x) * " \$\$")
Base.show(io::IO, ::MIME"text/latex", x::Symbolic) = print(io, "\$\$ " * latexify(x) * " \$\$")
Base.show(io::IO, ::MIME"text/latex", x::Equation) = print(io, "\$\$ " * latexify(x) * " \$\$")
Base.show(io::IO, ::MIME"text/latex", x::Vector{Equation}) = print(io, "\$\$ " * latexify(x) * " \$\$")
Base.show(io::IO, ::MIME"text/latex", x::AbstractArray{<:RCNum}) = print(io, "\$\$ " * latexify(x) * " \$\$")

_toexpr(O::ArrayOp) = _toexpr(O.term)

# `_toexpr` is only used for latexify
function _toexpr(O)
    if ismul(O)
        m = O
        numer = Any[]
        denom = Any[]

        # We need to iterate over each term in m, ignoring the numeric coefficient.
        # This iteration needs to be stable, so we can't iterate over m.dict.
        for term in Iterators.drop(sorted_arguments(m), isone(m.coeff) ? 0 : 1)
            if !ispow(term)
                push!(numer, _toexpr(term))
                continue
            end

            base = term.base
            pow  = term.exp
            isneg = (pow isa Number && pow < 0) || (iscall(pow) && operation(pow) === (-) && length(arguments(pow)) == 1)
            if !isneg
                if _isone(pow)
                    pushfirst!(numer, _toexpr(base))
                else
                    pushfirst!(numer, Expr(:call, :^, _toexpr(base), _toexpr(pow)))
                end
            else
                newpow = -1*pow
                if _isone(newpow)
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
    if issym(O) 
        sym = string(nameof(O))
        sym = replace(sym, NAMESPACE_SEPARATOR => ".")
        if length(sym) > 1
            sym = string("\\mathtt{", sym, "}")
        end
        return Symbol(sym)
    end
    !iscall(O) && return O

    op = operation(O)
    args = sorted_arguments(O)

    if (op===(*)) && (args[1] === -1)
        arg_mul = Expr(:call, :(*), _toexpr(args[2:end])...)
        return Expr(:call, :(-), arg_mul)
    end

    if op isa Differential
        num = args[1]
        den = op.x
        deg = 1
        while num isa Term && num.f isa Differential
            deg += 1
            den *= num.f.x
            num = first(arguments(num))
        end
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
        return Expr(:call, _toexpr(op), _toexpr(args)...)
    elseif op === getindex && symtype(args[1]) <: AbstractArray
        return getindex_to_symbol(O)
    elseif op === (\)
        return :(solve($(_toexpr(args[1])), $(_toexpr(args[2]))))
    elseif issym(op) && symtype(op) <: AbstractArray
        return :(_textbf($(nameof(op))))
    elseif op === identity
        return _toexpr(only(args)) # suppress identity transformations (e.g. "identity(π)" -> "π")
    end
    return Expr(:call, Symbol(op), _toexpr(args)...)
end
_toexpr(x::Integer) = x
_toexpr(x::AbstractFloat) = x

function _toexpr(eq::Equation)
    Expr(:(=), _toexpr(eq.lhs), _toexpr(eq.rhs))
end

_toexpr(eqs::AbstractArray) = map(eq->_toexpr(eq), eqs)
_toexpr(x::Num) = _toexpr(value(x))

function getindex_to_symbol(t)
    @assert iscall(t) && operation(t) === getindex && symtype(sorted_arguments(t)[1]) <: AbstractArray
    args = sorted_arguments(t)
    idxs = args[2:end]
    return :($(_toexpr(args[1]))[$(idxs...)])
end

function diffdenom(e)
    if issym(e)
        LaTeXString("\\mathrm{d}$e")
    elseif ispow(e)
        LaTeXString("\\mathrm{d}$(e.base)$(isone(e.exp) ? "" : "^{$(e.exp)}")")
    elseif ismul(e)
        LaTeXString(prod(
                "\\mathrm{d}$(k)$(isone(v) ? "" : "^{$v}")"
                for (k, v) in e.dict
               ))
    else
        e
    end
end
