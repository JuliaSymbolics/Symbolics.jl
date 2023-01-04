# define a context type
Cassette.@context SymbolicContext

# this is a near copy of `wrap_func_expr` in file "wrapper_types.jl"
# Instead of overloading `fname`, however, we specialize 
# `Cassette.overdub` on the symbolic argument types.
function ctx_wrap_func_expr(mod, expr)
    @assert expr.head == :function || (expr.head == :(=) &&
                                       expr.args[1] isa Expr &&
                                       expr.args[1].head == :call)

    def = splitdef(expr)

    sig = expr.args[1]
    body = def[:body]

    fname = def[:name]
    args = get(def, :args, [])
    kwargs = get(def, :kwargs, [])

    impl_name = Symbol(fname,"_", hash(string(args)*string(kwargs)))

    function kwargname(kwarg)
        if kwarg isa Expr && kwarg.head == :kw
            kwarg.args[1]
        elseif kwarg isa Expr && kwarg.head == :(...)
            kwarg.args[1]
        else
            kwarg
        end
    end

    function argname(arg)
        if arg isa Expr && (arg.head == :(::) || arg.head == :(...))
            arg.args[1]
        elseif arg isa Expr
            error("$arg not supported as an argument")
        else
            arg
        end
    end

    names = vcat(argname.(args), kwargname.(kwargs))

    function type_options(arg)
        if arg isa Expr && arg.head == :(::)
            T = Base.eval(mod, arg.args[2])
            has_symwrapper(T) ? (T, :(SymbolicUtils.Symbolic{<:$T}), wrapper_type(T)) :
                                (T,:(SymbolicUtils.Symbolic{<:$T}))
        elseif arg isa Expr && arg.head == :(...)
            Ts = type_options(arg.args[1])
            map(x->Vararg{x},Ts)
        else
            (Any,)
        end
    end

    types = map(type_options, args)

    impl = :(function $impl_name($(names...))
        println("dubbed")
        $body
    end)
    # TODO: maybe don't drop first lol
    methods = map(Iterators.drop(Iterators.product(types...), 1)) do Ts
        method_args = map(names, Ts) do n, T
            :($n::$T)
        end

        fbody = :(if any($iswrapped, ($(names...),))
                      $wrap($impl_name($([:($unwrap($arg)) for arg in names]...)))
                  else
                      $impl_name($(names...))
                  end)

        if isempty(kwargs)
            :(function Cassette.overdub(::SymbolicContext, ::typeof($(fname)), $(method_args...))
                  $fbody
              end)
        else
            :(function Cassette.overdub(::SymbolicContext, ::typeof($(fname)), $(method_args...); $(kwargs...))
                  $fbody
              end)
        end
    end

    quote
        $impl
        $(methods...)
    end |> esc
end

macro ctx_wrapped(expr)
    ctx_wrap_func_expr(__module__, expr)
end

# due to issue 
# https://github.com/JuliaPackaging/Requires.jl/issues/86
# we do all the wrapping in this file, too.

# "array-lib.jl"
@ctx_wrapped function getindex_posthook(f, x::AbstractArray)
    if hasmetadata(x, GetindexPosthookCtx)
        g = getmetadata(x, GetindexPosthookCtx)
        setmetadata(x,
                    GetindexPosthookCtx,
                    (res, args...) -> f(g(res, args...), args...))
    else
        setmetadata(x, GetindexPosthookCtx, f)
    end
end
@ctx_wrapped function -(x::CartesianIndex, y::CartesianIndex)
    CartesianIndex((tup(x) .- tup(y))...)
end
@ctx_wrapped function +(x::CartesianIndex, y::CartesianIndex)
    CartesianIndex((tup(x) .+ tup(y))...)
end
@ctx_wrapped function Base.adjoint(A::AbstractMatrix)
    @syms i::Int j::Int
    @arrayop (i, j) A[j, i] term=A'
end
@ctx_wrapped function Base.adjoint(b::AbstractVector)
    @syms i::Int
    @arrayop (1, i) b[i] term=b'
end
@ctx_wrapped (*)(A::AbstractMatrix, B::AbstractMatrix) = _matmul(A, B)
@ctx_wrapped (*)(A::AbstractVector, B::AbstractMatrix) = _matmul(A, B)
@ctx_wrapped (*)(A::AbstractMatrix, b::AbstractVector) = _matvec(A, b)
@ctx_wrapped Base.map(f, x::AbstractArray) = _map(f, x)
@ctx_wrapped Base.map(f, x::AbstractArray, xs...) = _map(f, x, xs...)
@ctx_wrapped Base.map(f, x, y::AbstractArray, z...) = _map(f, x, y, z...)
@ctx_wrapped Base.map(f, x, y, z::AbstractArray, w...) = _map(f, x, y, z, w...)
@ctx_wrapped function Base.mapreduce(f, g, x::AbstractArray; dims=:, kw...)
    idx = makesubscripts(ndims(x))
    out_idx = [dims == (:) || i in dims ? 1 : idx[i] for i = 1:ndims(x)]
    expr = f(x[idx...])
    T = symtype(g(expr, expr))
    if dims === (:)
        return Term{T}(_mapreduce, [f, g, x, dims, (kw...,)])
    end

    Atype = propagate_atype(_mapreduce, f, g, x, dims, (kw...,))
    ArrayOp(Atype{T, ndims(x)},
            (out_idx...,),
            expr,
            g,
            Term{Any}(_mapreduce, [f, g, x, dims, (kw...,)]))
end
for (ff, opts) in [sum => (identity, +, false),
                  prod => (identity, *, true),
                  any => (identity, (|), false),
                  all => (identity, (&), true)]

    f, g, init = opts
    @eval @ctx_wrapped function $(ff)(x::AbstractArray;
                                     dims=:, init=$init)
        mapreduce($f, $g, x, dims=dims, init=init)
    end
    @eval @ctx_wrapped function $(ff)(f::Function, x::AbstractArray;
                                     dims=:, init=$init)
        mapreduce(f, $g, x, dims=dims, init=init)
    end
end

# "arrays.jl"
@ctx_wrapped function Base.:(\)(A::AbstractMatrix, b::AbstractVecOrMat)
    t = arrterm(\, A, b)
    setmetadata(t, ScalarizeCache, Ref{Any}(nothing))
end
@ctx_wrapped function Base.inv(A::AbstractMatrix)
    t = arrterm(inv, A)
    setmetadata(t, ScalarizeCache, Ref{Any}(nothing))
end
@ctx_wrapped function LinearAlgebra.det(x::AbstractMatrix; laplace=true)
    Term{eltype(x)}(_det, [x, laplace])
end
@ctx_wrapped Base.isempty(x::AbstractArray) = shape(unwrap(x)) !== Unknown() && _iszero(length(x))