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
function Base.adjoint(A::SymMat)
    @syms i::Int j::Int
    @arrayop A' (i, j) A[j, i]
end

function Base.adjoint(b::SymVec)
    @syms i::Int
    @arrayop b' (1, i) b[i]
end

import Base: *, \

function isdot(A::Term, ::SymVec)
    operation(A) === (adjoint) && symtype(arguments(A)[1]) <: AbstractVector
end
isdot(A::ArrayOp, b::SymVec) = isdot(A.term, b)
isdot(A, b) = false

function (*)(A::SymMat, B::SymMat)
    @syms i::Int j::Int k::Int
    @arrayop (A*B) (i, j) A[i, k] * B[k, j]
end

function (*)(A::SymMat, b::SymVec)
    if isdot(A, b)
        return Term{Number}(*, [A, b])
    end
    @syms i::Int k::Int
    @arrayop (A*b) (i,) A[i, k] * b[k]
end

#################### MAP-REDUCE ################
#

Base.map(f, x::SymArray) = _map(f, x)
Base.map(f, x::SymArray, xs...) = _map(f, x, xs...)
Base.map(f, x, y::SymArray, z...) = _map(f, x, y, z...)
Base.map(f, x, y, z::SymArray, w...) = _map(f, x, y, z, w...)

function _map(f, x, xs...)
    N = ndims(x)
    idx = makesubscripts(N)

    expr = f(map(a->a[idx...], [x, xs...])...)

    Atype = propagate_atype(map, f, x, xs...)
    ArrayOp(Atype{symtype(expr), N},
            (idx...,),
            f(expr),
            +,
            Term{Any}(map, [f, x, xs...]))
end

@inline _mapreduce(f, x, dims, kw) = mapreduce(f, x; dims=dims, kw...)

function Base.mapreduce(f, g, x::SymArray; dims=:, kw...)
    if dims === (:)
        return Term{Number}(_mapreduce, [f, g, x, dims, (kw...,)])
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
    @eval function (::$(typeof(ff)))(x::SymArray;
                                     dims=:, init=$init)
        mapreduce($f, $g, x, dims=dims, init=init)
    end
end
