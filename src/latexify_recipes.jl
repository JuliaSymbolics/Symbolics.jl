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
    rhs = prettify_expr.(_toexpr(rhs; canonicalize=false))
    rhs = [postwalk(x -> x isa Expr && length(x.args) == 1 ? x.args[1] : x, eq) for eq in rhs]
    rhs = [postwalk(x -> x isa Expr && x.args[1] == :_derivative && length(x.args[2].args) == 2 ? :($(Symbol(:d, x.args[2]))/($(Symbol(:d, x.args[2].args[2])))) : x, eq) for eq in rhs]
    rhs = [postwalk(x -> x isa Expr && x.args[1] == :_derivative ? "\\frac{d\\left($(Latexify.latexraw(x.args[2]))\\right)}{d$(Latexify.latexraw(x.args[3]))}" : x, eq) for eq in rhs]

    lhs = getfield.(eqs, :lhs)
    lhs = prettify_expr.(_toexpr(lhs; canonicalize=false))
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
    return _toexpr(n, canonicalize=false)
end

# `_toexpr` is only used for latexify
function _toexpr(O; canonicalize=true)
    if canonicalize
        canonical, O = canonicalexpr(O)
        canonical && return O
    else
        !istree(O) && return O
    end
    op = operation(O)
    args = arguments(O)
    
    if (op===(*)) && (args[1] === -1)
        println("here")
	println(args)
    	arg_mul = Expr(:call, :(*), _toexpr(args[2:end]; canonicalize=canonicalize)...)
        return Expr(:call, :(-), arg_mul)
    end
    
    if op isa Differential
        ex = _toexpr(args[1]; canonicalize=canonicalize)
        wrt = _toexpr(op.x; canonicalize=canonicalize)
        return :(_derivative($ex, $wrt))
    elseif op isa Sym
        isempty(args) && return nameof(op)
        return Expr(:call, _toexpr(op; canonicalize=canonicalize), _toexpr(args; canonicalize=canonicalize)...)
    end
    return Expr(:call, Symbol(op), _toexpr(args; canonicalize=canonicalize)...)
end
_toexpr(s::Sym; kw...) = nameof(s)

function canonicalexpr(O)
    !istree(O) && return true, O
    op = operation(O)
    args = arguments(O)
    if op === (^)
        if length(args) == 2 && args[2] isa Number && args[2] < 0
            ex = _toexpr(args[1])
            if args[2] == -1
                expr = Expr(:call, :inv, ex)
            else
                expr = Expr(:call, :^, Expr(:call, inv, ex), -args[2])
            end
            return true, expr
        end
    end
    return false, O
end
for fun in [:_toexpr]
    @eval begin
        function $fun(eq::Equation; kw...)
            Expr(:(=), $fun(eq.lhs; kw...), $fun(eq.rhs; kw...))
        end

        $fun(eqs::AbstractArray; kw...) = map(eq->$fun(eq; kw...), eqs)
        $fun(x::Integer; kw...) = x
        $fun(x::AbstractFloat; kw...) = x
    end
end
_toexpr(x::Num; kw...) = _toexpr(value(x); kw...)
