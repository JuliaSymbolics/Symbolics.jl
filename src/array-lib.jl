##### getindex #####
function Base.getindex(x::ArrayLike, idx...)
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

makesubscripts(n) = [Sym{Int}(Symbol("i_$i")) for i in 1:n]

using Base.Broadcast

Base.broadcastable(s::ArrayLike) = s
struct SymBroadcast <: Broadcast.BroadcastStyle end
Broadcast.BroadcastStyle(::Type{<:ArrayLike}) = SymBroadcast()
Broadcast.result_style(::SymBroadcast) = SymBroadcast()
Broadcast.BroadcastStyle(::SymBroadcast, ::Broadcast.BroadcastStyle) = SymBroadcast()

isonedim(x, i) = shape(x) == Unknown() ? false : isone(size(x, i))

function Broadcast.materialize(bc::Broadcast.Broadcasted{SymBroadcast})
    # Do the thing here
    ndim = mapfoldl(ndims, max, bc.args, init=0)
    subscripts = makesubscripts(ndim)

    expr_args′ = map(enumerate(bc.args)) do (i, x)
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
    ArrayOp{Atype{symtype(expr), ndim}}(
        (subscripts...,),
        expr,
        +,
        Term{Any}(broadcast, [bc.f, bc.args...]))
end

#################### TRANSPOSE ################
#
function Base.adjoint(A::ArrayLike)
    @syms i::Int j::Int

    N = ndims(A)
    if N == 1
        @arrayop A' (1, j) A[j]
    elseif N == 2
        @arrayop A' (i, j) A[j, i]
    else
        error("Can adjoint only a vector or a matrix")
    end
end
