import Base: getindex

##### getindex #####
function Base.getindex(x::SymArray, idx...)
    if all(i->symtype(i) <: Integer, idx)
        Term{eltype(symtype(x))}(getindex, [x, idx...])
    else
        arrterm(getindex, x, idx...)
    end
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


#### Broadcast ####
#

function makesubscripts(n)
    set = 'i':'z'
    m = length(set)
    map(1:n) do i
        repeats = ceil(Int, i / m)
        c = set[(i-1) % m + 1]
        Sym{Int}(Symbol(join([c for _ in 1:repeats], "")))
    end
end

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

    expr_args′ = map(bc.args) do x
        if ndims(x) != 0
            subs = map(i-> isonedim(x, i) ?
                       1 : subscripts[i], 1:ndims(x))
            x[subs...]
        else
            x
        end
    end

    expr = term(bc.f, expr_args′...)
    Atype = propagate_atype(broadcast, bc.f, bc.args...)
    ArrayOp(Atype{symtype(expr), ndim},
            (subscripts...,),
            expr,
            +,
            Term{Any}(broadcast, [bc.f, bc.args...]))
end

#################### TRANSPOSE ################
#
@wrapped function Base.adjoint(A::AbstractMatrix)
    @syms i::Int j::Int
    @arrayop A' (i, j) A[j, i]
end

@wrapped function Base.adjoint(b::AbstractVector)
    @syms i::Int
    @arrayop b' (1, i) b[i]
end

import Base: *, \

using LinearAlgebra

isdot(A::Adjoint, ::SymVec) = true

# TODO: add more such methods
function getindex(A::AbstractArray, i::Symbolic{<:Integer}...)
    Term{eltype(A)}(getindex, [A, i...])
end

function getindex(A::AbstractArray, i::Int, j::Symbolic{<:Integer})
    Term{eltype(A)}(getindex, [A, i, j])
end

function getindex(A::AbstractArray, j::Symbolic{<:Integer}, i::Int)
    Term{eltype(A)}(getindex, [A, j, i])
end

function isdot(A::Term, ::SymVec)
    operation(A) === (adjoint) && symtype(arguments(A)[1]) <: AbstractVector
end
isdot(A::ArrayOp, b::Union{SymVec, Vector}) = isdot(A.term, b)
isdot(A, b) = false

function _matmul(A,B)
    @syms i::Int j::Int k::Int
    @arrayop (A*B) (i, j) A[i, k] * B[k, j]
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
    @arrayop (A*b) (i,) A[i, k] * b[k]
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

@inline _mapreduce(f, x, dims, kw) = mapreduce(f, x; dims=dims, kw...)

@wrapped function Base.mapreduce(f, g, x::AbstractArray; dims=:, kw...)
    if dims === (:)
        T = SymbolicUtils.promote_symtype(f, eltype(x))
        S = SymbolicUtils.promote_symtype(g, T, T)
        return Term{S}(_mapreduce, [f, g, x, dims, (kw...,)])
    end

    idx = makesubscripts(ndims(x))
    out_idx = [i in dims ? 1 : idx[i] for i = 1:ndims(x)]
    expr = f(x[idx...])

    Atype = propagate_atype(_mapreduce, f, g, x, dims, (kw...,))
    ArrayOp(Atype{symtype(expr), ndims(x)},
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
end


# Wrapped array should wrap the elements too
function Base.getindex(x::Arr, idx...)
    wrap(unwrap(x)[idx...])
end
function Base.getindex(x::Arr, idx::Symbolic{<:Integer}...)
    wrap(unwrap(x)[idx...])
end

struct SymWrapBroadcast <: Broadcast.BroadcastStyle end

Base.broadcastable(s::Arr) = s

Broadcast.BroadcastStyle(::Type{<:Arr}) = SymWrapBroadcast()

Broadcast.result_style(::SymWrapBroadcast) = SymWrapBroadcast()

Broadcast.BroadcastStyle(::SymWrapBroadcast,
                         ::Broadcast.BroadcastStyle) = SymWrapBroadcast()

function Broadcast.materialize(bc::Broadcast.Broadcasted{SymWrapBroadcast})
    wrap(broadcast(bc.f, unwrap.(bc.args)...))
end

