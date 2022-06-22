import Base: getindex

##### getindex #####
struct GetindexPosthookCtx end

@wrapped function getindex_posthook(f, x::AbstractArray)
    if hasmetadata(x, GetindexPosthookCtx)
        g = getmetadata(x, GetindexPosthookCtx)
        setmetadata(x,
                    GetindexPosthookCtx,
                    (res, args...) -> f(g(res, args...), args...))
    else
        setmetadata(x, GetindexPosthookCtx, f)
    end
end
function Base.getindex(x::SymArray, idx::CartesianIndex)
    return x[Tuple(idx)...]
end

function Base.getindex(x::SymArray, idx...)
    idx = unwrap.(idx)
    if shape(x) !== Unknown() && all(i->i isa Integer, idx)
        II = CartesianIndices(axes(x))
        @boundscheck begin
            if !checkbounds(Bool, II, idx...)
                throw(BoundsError(x, idx))
            end
        end
        ii = II[idx...]
        res = Term{eltype(symtype(x))}(getindex, [x, Tuple(ii)...])
    elseif all(i->symtype(i) <: Integer, idx)
        shape(x) !== Unknown() && @boundscheck begin
            if length(idx) > 1
                for (a, i) in zip(axes(x), idx)
                    if i isa Integer && !(i in a)
                        throw(BoundsError(x, idx))
                    end
                end
            end
        end
        res = Term{eltype(symtype(x))}(getindex, [x, idx...])
    elseif length(idx) == 1 && symtype(first(idx)) <: CartesianIndex
        i = first(idx)
        ii = i isa CartesianIndex ? Tuple(i) : arguments(i)

        return getindex(x, ii...)
    else
        input_idx = []
        output_idx = []
        ranges = Dict{Sym, AbstractRange}()
        subscripts = makesubscripts(length(idx))
        for (j, i) in enumerate(idx)
            if symtype(i) <: Integer
                push!(input_idx, i)
            elseif i isa Colon
                push!(output_idx, subscripts[j])
                push!(input_idx, subscripts[j])
            elseif i isa AbstractVector
                isym = subscripts[j]
                push!(output_idx, isym)
                push!(input_idx, isym)
                ranges[isym] = i
            else
                error("Don't know how to index by $i")
            end
        end

        term = Term{Any}(getindex, [x, idx...])
        T = eltype(symtype(x))
        N = ndims(x) - count(i->symtype(i) <: Integer, idx)
        res = ArrayOp(atype(symtype(x)){T,N},
                      (output_idx...,),
                      x[input_idx...],
                      +,
                      term,
                      ranges)
    end

    if hasmetadata(x, GetindexPosthookCtx)
        f = getmetadata(x, GetindexPosthookCtx, identity)
        f(res, x, idx...)
    else
        res
    end
end

# Wrapped array should wrap the elements too
function Base.getindex(x::Arr, idx...)
    wrap(unwrap(x)[idx...])
end
function Base.getindex(x::Arr, idx::Symbolic{<:Integer}...)
    wrap(unwrap(x)[idx...])
end

function Base.CartesianIndex(x::Symbolic{<:Integer}, xs::Symbolic{<:Integer}...)
    term(CartesianIndex, x, xs..., type=CartesianIndex)
end

import Base: +, -
tup(c::CartesianIndex) = Tuple(c)
tup(c::Term{CartesianIndex}) = arguments(c)
@wrapped function -(x::CartesianIndex, y::CartesianIndex)
    CartesianIndex((tup(x) .- tup(y))...)
end
@wrapped function +(x::CartesianIndex, y::CartesianIndex)
    CartesianIndex((tup(x) .+ tup(y))...)
end

function propagate_ndims(::typeof(getindex), x, idx...)
    ndims(x) - count(x->symtype(x) <: Integer, idx)
end

function propagate_shape(::typeof(getindex), x, idx...)
    @oops axes = shape(x)

    idx1 = to_indices(CartesianIndices(axes), axes, idx)
    ([1:length(x) for x in idx1 if !(symtype(x) <: Number)]...,)
end

propagate_eltype(::typeof(getindex), x, idx...) = geteltype(x)

function SymbolicUtils.promote_symtype(::typeof(getindex), X, ii...)
    @assert all(i-> i <:Integer, ii) "user arrterm to create arr term."

    eltype(X)
end


#### Broadcast ####
#

using Base.Broadcast

Base.broadcastable(s::SymArray) = s
struct SymBroadcast <: Broadcast.BroadcastStyle end
Broadcast.BroadcastStyle(::Type{<:SymArray}) = SymBroadcast()
Broadcast.result_style(::SymBroadcast) = SymBroadcast()
Broadcast.BroadcastStyle(::SymBroadcast, ::Broadcast.BroadcastStyle) = SymBroadcast()

isonedim(x, i) = shape(x) == Unknown() ? false : isone(size(x, i))

function Broadcast.materialize(bc::Broadcast.Broadcasted{SymBroadcast})
    # Do the thing here
    ndim = mapfoldl(ndims, max, bc.args, init=0)
    subscripts = makesubscripts(ndim)

    onedim_count = mapreduce(+, bc.args) do x
        if ndims(x) != 0
            map(i-> isonedim(x, i) ? 1 : 0, 1:ndim)
        else
            map(i-> 1, 1:ndim)
        end
    end

    extruded = map(x->x < length(bc.args), onedim_count)

    expr_args′ = map(bc.args) do x
        if ndims(x) != 0
            subs = map(i-> extruded[i] && isonedim(x, i) ?
                       1 : subscripts[i], 1:ndims(x))
            x[subs...]
        else
            x
        end
    end

    expr = term(bc.f, expr_args′...) # Imagine x .=> y -- if you don't have a term
                                     # then you get pairs, and index matcher cannot
                                     # recurse into pairs
    Atype = propagate_atype(broadcast, bc.f, bc.args...)
    ArrayOp(Atype{symtype(expr), ndim},
            (subscripts...,),
            expr,
            +,
            Term{Any}(broadcast, [bc.f, bc.args...]))
end

# On wrapper:
struct SymWrapBroadcast <: Broadcast.BroadcastStyle end

Base.broadcastable(s::Arr) = s

Broadcast.BroadcastStyle(::Type{<:Arr}) = SymWrapBroadcast()

Broadcast.result_style(::SymWrapBroadcast) = SymWrapBroadcast()

Broadcast.BroadcastStyle(::SymWrapBroadcast,
                         ::Broadcast.BroadcastStyle) = SymWrapBroadcast()
Broadcast.BroadcastStyle(::SymBroadcast,
                         ::SymWrapBroadcast) = Broadcast.Unknown()

function Broadcast.materialize(bc::Broadcast.Broadcasted{SymWrapBroadcast})
    args = map(bc.args) do arg
        if arg isa Broadcast.Broadcasted
            return Broadcast.materialize(arg)
        else
            return arg
        end
    end
    wrap(broadcast(bc.f, map(unwrap, args)...))
end


#################### TRANSPOSE ################
#
@wrapped function Base.adjoint(A::AbstractMatrix)
    @syms i::Int j::Int
    @arrayop (i, j) A[j, i] term=A'
end

@wrapped function Base.adjoint(b::AbstractVector)
    @syms i::Int
    @arrayop (1, i) b[i] term=b'
end

import Base: *, \

using LinearAlgebra

isdot(A, b) = isadjointvec(A) && ndims(b) == 1

isadjointvec(A::Adjoint) = ndims(parent(A)) == 1
isadjointvec(A::Transpose) = ndims(parent(A)) == 1

function isadjointvec(A::Term)
    (operation(A) === (adjoint) ||
     operation(A) == (transpose)) && ndims(arguments(A)[1]) == 1
end

isadjointvec(A::ArrayOp) = isadjointvec(A.term)

isadjointvec(A) = false

# TODO: add more such methods
function getindex(A::AbstractArray, i::Symbolic{<:Integer}, ii::Symbolic{<:Integer}...)
    Term{eltype(A)}(getindex, [A, i, ii...])
end

function getindex(A::AbstractArray, i::Int, j::Symbolic{<:Integer})
    Term{eltype(A)}(getindex, [A, i, j])
end

function getindex(A::AbstractArray, j::Symbolic{<:Integer}, i::Int)
    Term{eltype(A)}(getindex, [A, j, i])
end

function getindex(A::Arr, i::Int, j::Symbolic{<:Integer})
    wrap(unwrap(A)[i,j])
end

function getindex(A::Arr, j::Symbolic{<:Integer}, i::Int)
    wrap(unwrap(A)[j,i])
end

function _matmul(A,B)
    @syms i::Int j::Int k::Int
    if isadjointvec(A)
        op = operation(A.term)
        return op(op(B) * first(arguments(A.term)))
    end
    return @arrayop (i, j) A[i, k] * B[k, j] term=(A*B)
end

@wrapped (*)(A::AbstractMatrix, B::AbstractMatrix) = _matmul(A, B)
@wrapped (*)(A::AbstractVector, B::AbstractMatrix) = _matmul(A, B)

function _matvec(A,b)
    @syms i::Int k::Int
    if isdot(A, b)
        make_shape((i,), A[i, k] * b[k]) # This is a dimension check
        T = SymbolicUtils.promote_symtype(*, eltype(A), eltype(b))
        S = SymbolicUtils.promote_symtype(+, T,T)
        return Term{S}(*, [A, b])
    end
    @arrayop (i,) A[i, k] * b[k] term=(A*b)
end
@wrapped (*)(A::AbstractMatrix, b::AbstractVector) = _matvec(A, b)
(*)(A::AbstractArray, b::Union{AbstractVector, SymVec}) = _matvec(A, b)

#################### MAP-REDUCE ################
#

@wrapped Base.map(f, x::AbstractArray) = _map(f, x)
@wrapped Base.map(f, x::AbstractArray, xs...) = _map(f, x, xs...)
@wrapped Base.map(f, x, y::AbstractArray, z...) = _map(f, x, y, z...)
@wrapped Base.map(f, x, y, z::AbstractArray, w...) = _map(f, x, y, z, w...)

function _map(f, x, xs...)
    N = ndims(x)
    idx = makesubscripts(N)

    expr = f(map(a->a[idx...], [x, xs...])...)

    Atype = propagate_atype(map, f, x, xs...)
    ArrayOp(Atype{symtype(expr), N},
            (idx...,),
            expr,
            +,
            Term{Any}(map, [f, x, xs...]))
end

@inline _mapreduce(f, g, x, dims, kw) = mapreduce(f, g, x; dims=dims, kw...)

function scalarize_op(::typeof(_mapreduce), t)
    f,g,x,dims,kw = arguments(t)
    # we wrap and unwrap to make things work smoothly.
    # we need the result unwrapped to allow recursive scalarize to work.
    unwrap(_mapreduce(f, g, collect(wrap(x)), dims, kw))
end

@wrapped function Base.mapreduce(f, g, x::AbstractArray; dims=:, kw...)
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
    @eval @wrapped function (::$(typeof(ff)))(x::AbstractArray;
                                     dims=:, init=$init)
        mapreduce($f, $g, x, dims=dims, init=init)
    end
    @eval @wrapped function (::$(typeof(ff)))(f::Function, x::AbstractArray;
                                     dims=:, init=$init)
        mapreduce(f, $g, x, dims=dims, init=init)
    end
end
