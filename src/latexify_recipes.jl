prettify_expr(expr) = expr
prettify_expr(f::Function) = nameof(f)
prettify_expr(expr::Expr) = Expr(expr.head, prettify_expr.(expr.args)...)

@latexrecipe function f(eqs::Vector{Equation})
    # Set default option values.
    env --> :align
    cdot --> false

    # Convert both the left- and right-hand sides to expressions of basic types
    # that latexify can deal with

    rhs = getfield.(eqs, :rhs)
    rhs = prettify_expr.(_toexpr(rhs))
    rhs = [postwalk(x -> x isa Expr && length(x.args) == 1 ? x.args[1] : x, eq) for eq in rhs]
    rhs = [postwalk(x -> x isa Expr && x.args[1] == :_derivative && length(x.args[2].args) == 2 ? :($(Symbol(:d, x.args[2]))/($(Symbol(:d, x.args[2].args[2])))) : x, eq) for eq in rhs]
    rhs = [postwalk(x -> x isa Expr && x.args[1] == :_derivative ? "\\frac{d\\left($(Latexify.latexraw(x.args[2]))\\right)}{d$(Latexify.latexraw(x.args[3]))}" : x, eq) for eq in rhs]

    lhs = getfield.(eqs, :lhs)
    lhs = prettify_expr.(_toexpr(lhs))
    lhs = [postwalk(x -> x isa Expr && length(x.args) == 1 ? x.args[1] : x, eq) for eq in lhs]
    lhs = [postwalk(x -> x isa Expr && x.args[1] == :_derivative && length(x.args[2].args) == 2 ? :($(Symbol(:d, x.args[2]))/($(Symbol(:d, x.args[2].args[2])))) : x, eq) for eq in lhs]
    lhs = [postwalk(x -> x isa Expr && x.args[1] == :_derivative ? "\\frac{d\\left($(Latexify.latexraw(x.args[2]))\\right)}{d$(Latexify.latexraw(x.args[3]))}" : x, eq) for eq in lhs]

    return lhs, rhs
end

Base.show(io::IO, ::MIME"text/latex", x::Num) = print(io, latexify(x))
Base.show(io::IO, ::MIME"text/latex", x::Symbolic) = print(io, latexify(x))
Base.show(io::IO, ::MIME"text/latex", x::Vector{Equation}) = print(io, latexify(x))
Base.show(io::IO, ::MIME"text/latex", x::AbstractArray{Num}) = print(io, latexify(x))

@latexrecipe function f(n::Num)
    return _toexpr(n)
end

# `_toexpr` is only used for latexify
function _toexpr(O)
    !istree(O) && return O

    op = operation(O)
    args = arguments(O)

    if (op===(*)) && (args[1] === -1)
        arg_mul = Expr(:call, :(*), _toexpr(args[2:end])...)
        return Expr(:call, :(-), arg_mul)
    end

    if op isa Differential
        ex = _toexpr(args[1])
        wrt = _toexpr(op.x)
        return :(_derivative($ex, $wrt))
    elseif op isa Sym
        isempty(args) && return nameof(op)
        return Expr(:call, _toexpr(op), _toexpr(args)...)
    end
    return Expr(:call, Symbol(op), _toexpr(args)...)
end
function _toexpr(m::Mul{<:Number})
    numer = Any[]
    denom = Any[]

    for (base, pow) in m.dict
        if pow > 0
            if isone(pow)
                pushfirst!(numer, _toexpr(base))
            else
                pushfirst!(numer, Expr(:call, :^, _toexpr(base), pow))
            end
        else
            if isone(-1*pow)
                pushfirst!(denom, _toexpr(base))
            else
                pushfirst!(denom, Expr(:call, :^, _toexpr(base), -1*pow))
            end
        end
    end

    if isempty(numer) || !isone(abs(m.coeff))
        numer_expr = Expr(:call, :*, abs(m.coeff), numer...)
    else
        numer_expr = Expr(:call, :*, numer...)
    end

    if isempty(denom)
        frac_expr = numer_expr
    else
        denom_expr = Expr(:call, :*, denom...)
        frac_expr = Expr(:call, :/, numer_expr, denom_expr)
    end

    if m.coeff < 0
        return Expr(:call, :-, frac_expr)
    else
        return frac_expr
    end
end
_toexpr(s::Sym) = nameof(s)
_toexpr(x::Integer) = x
_toexpr(x::AbstractFloat) = x

function _toexpr(eq::Equation)
    Expr(:(=), _toexpr(eq.lhs), _toexpr(eq.rhs))
end

_toexpr(eqs::AbstractArray) = map(eq->_toexpr(eq), eqs)
_toexpr(x::Num) = _toexpr(value(x))
