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

recipe(n; kw...) = latexify_derivatives(cleanup_exprs(_toexpr(n; kw...)))

@latexrecipe function f(n::Num)
    # SMALL HACK: forward any keyword arguments passed to latexify()
    # must be done *before* attributes are set (since they modify kwargs internally inside Latexify)
    # (see https://github.com/korsbo/Latexify.jl/issues/320)
    ret = recipe(n; kwargs...)

    env --> :equation
    cdot --> false
    fmt --> FancyNumberFormatter(5)
    index --> :subscript
    snakecase --> true
    safescripts --> true

    return ret
end

@latexrecipe function f(z::Complex{Num})
    # same hack as above
    if iszero(z.im)
        ret = :($(recipe(z.re; kwargs...)))
    elseif iszero(z.re)
        ret = :($(recipe(z.im; kwargs...)) * $im)
    else
        ret = :($(recipe(z.re; kwargs...)) + $(recipe(z.im; kwargs...)) * $im)
    end

    env --> :equation
    cdot --> false
    index --> :subscript

    return ret
end

@latexrecipe function f(n::ArrayOp)
    ret = recipe(n.term) # same hack as above

    env --> :equation
    cdot --> false
    index --> :subscript
    
    return ret
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
    ret = recipe(n; kwargs...) # same hack as above

    env --> :equation
    cdot --> false
    index --> :subscript

    return ret
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

    if eq.lhs isa Connection
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

# `_toexpr` is only used for latexify
function _toexpr(O; kw...)
    if ismul(O)
        m = O
        numer = Any[]
        denom = Any[]

        # We need to iterate over each term in m, ignoring the numeric coefficient.
        # This iteration needs to be stable, so we can't iterate over m.dict.
        for term in Iterators.drop(sorted_arguments(m), isone(m.coeff) ? 0 : 1)
            if !ispow(term)
                push!(numer, _toexpr(term; kw...))
                continue
            end

            base = term.base
            pow  = term.exp
            isneg = (pow isa Number && pow < 0) || (iscall(pow) && operation(pow) === (-) && length(arguments(pow)) == 1)
            if !isneg
                if _isone(pow)
                    pushfirst!(numer, _toexpr(base; kw...))
                else
                    pushfirst!(numer, Expr(:call, :^, _toexpr(base; kw...), _toexpr(pow; kw...)))
                end
            else
                newpow = -1*pow
                if _isone(newpow)
                    pushfirst!(denom, _toexpr(base; kw...))
                else
                    pushfirst!(denom, Expr(:call, :^, _toexpr(base; kw...), _toexpr(newpow; kw...)))
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

        if isreal(m.coeff) && m.coeff < 0
            return Expr(:call, :-, frac_expr)
        else
            return frac_expr
        end
    end
    if issym(O)
        name = String(nameof(O))
        name = replace(name, NAMESPACE_SEPARATOR => ".")
        return Symbol(_toexpr(name; kw...))
    end
    !iscall(O) && return O

    op = operation(O)
    args = sorted_arguments(O)

    if (op===(*)) && (args[1] === -1)
        arg_mul = Expr(:call, :(*), _toexpr(args[2:end]; kw...)...)
        return Expr(:call, :(-), arg_mul)
    end

    if op isa Differential
        num = args[1]
        den = op.x
        deg = 1
        while num isa Term && num.f isa Differential
            deg += 1
            den *= num.f.x
            num = num.arguments[1]
        end
        return :(_derivative($(_toexpr(num; kw...)), $den, $deg))
    elseif op isa Integral
        lower = op.domain.domain.left
        upper = op.domain.domain.right
        vars = op.domain.variables
        integrand = args[1]
        var = if vars isa Tuple
            Expr(:call, :(*), _toexpr(vars...; kw...))
        else
            _toexpr(vars; kw...)
        end
        return Expr(:call, :_integral, _toexpr(lower; kw...), _toexpr(upper; kw...), vars, _toexpr(integrand; kw...))
    elseif symtype(op) <: FnType
        isempty(args) && return nameof(op)
        return Expr(:call, _toexpr(op; kw...), _toexpr(args; kw...)...)
    elseif op === getindex && symtype(args[1]) <: AbstractArray
        return getindex_to_symbol(O)
    elseif op === (\)
        return :(solve($(_toexpr(args[1]; kw...)), $(_toexpr(args[2]; kw...))))
    elseif issym(op) && symtype(op) <: AbstractArray
        return :(_textbf($(nameof(op))))
    elseif op === identity
        return _toexpr(only(args); kw...) # suppress identity transformations (e.g. "identity(π)" -> "π")
    end
    return Expr(:call, Symbol(op), _toexpr(args; kw...)...)
end
_toexpr(x::Integer; kw...) = x
_toexpr(x::AbstractFloat; kw...) = x
_toexpr(x::AbstractString; variable_formatter = identity) = variable_formatter(x) # variable_formatter can modify latexified variable names
_toexpr(eq::Equation; kw...) = Expr(:(=), _toexpr(eq.lhs; kw...), _toexpr(eq.rhs; kw...))
_toexpr(eqs::AbstractArray; kw...) = map(eq->_toexpr(eq; kw...), eqs)
_toexpr(O::ArrayOp; kw...) = _toexpr(O.term; kw...)
_toexpr(x::Num; kw...) = _toexpr(value(x); kw...)

function getindex_to_symbol(t)
    @assert iscall(t) && operation(t) === getindex && symtype(sorted_arguments(t)[1]) <: AbstractArray
    args = sorted_arguments(t)
    idxs = args[2:end]
    return :($(_toexpr(args[1]; kw...))[$(idxs...)])
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
