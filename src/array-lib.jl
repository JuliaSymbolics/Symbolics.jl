import Base: getindex
inner_unwrap(x) = x isa AbstractArray ? unwrap.(x) : x

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

# Wrapped array should wrap the elements too
function Base.getindex(x::Arr, idx...)
    wrap(unwrap(x)[idx...])
end
const SymIdxT = Union{Num, BasicSymbolic{VartypeT}}
function Base.getindex(x::Arr, idx::SymIdxT, idxs...)
    wrap(unwrap(x)[idx, idxs...])
end
function Base.getindex(x::Arr, i1, idx::SymIdxT, idxs...)
    wrap(unwrap(x)[i1, idx, idxs...])
end
function Base.getindex(x::Arr, i1::SymIdxT, idx::SymIdxT, idxs...)
    wrap(unwrap(x)[i1, idx, idxs...])
end
function Base.getindex(x::Arr, i1, i2, idx::SymIdxT, idxs...)
    wrap(unwrap(x)[i1, i2, idx, idxs...])
end
function Base.getindex(x::Arr, i1, i2::SymIdxT, idx::SymIdxT, idxs...)
    wrap(unwrap(x)[i1, i2, idx, idxs...])
end
function Base.getindex(x::Arr, i1::SymIdxT, i2::SymIdxT, idx::SymIdxT, idxs...)
    wrap(unwrap(x)[i1, i2, idx, idxs...])
end

import Base: +, -, *

#### Broadcast ####

# On wrapper:
struct SymWrapBroadcast <: Broadcast.BroadcastStyle end

Broadcast.broadcastable(s::Arr) = s

Broadcast.BroadcastStyle(::Type{<:Arr}) = SymWrapBroadcast()

Broadcast.result_style(::SymWrapBroadcast) = SymWrapBroadcast()

Broadcast.BroadcastStyle(::SymWrapBroadcast,
    ::Broadcast.BroadcastStyle) = SymWrapBroadcast()
Broadcast.BroadcastStyle(::SymbolicUtils.SymBroadcast,
    ::SymWrapBroadcast) = Broadcast.Unknown()

unwrap_broadcasts(head, args...) = (unwrap_broadcast(head), unwrap_broadcasts(args...)...)
unwrap_broadcasts() = ()
unwrap_broadcast(x) = unwrap(x)
function unwrap_broadcast(bc::Broadcast.Broadcasted{SymWrapBroadcast})
    Broadcast.Broadcasted{SymbolicUtils.SymBroadcast{VartypeT}}(bc.f, unwrap_broadcasts(bc.args...), bc.axes)
end

function Broadcast.copy(bc::Broadcast.Broadcasted{SymWrapBroadcast})
    return wrap(copy(unwrap_broadcast(bc)))
end

#################### POLYADIC ################

function *(x::Arr, args...)
    return wrap(*(unwrap(x), args...))
end

function *(a::SymbolicUtils.PolyadicNumericOpFirstArgT, b::Arr, bs...)
    return wrap(*(a, unwrap(b), bs...))
end
function *(a::LinearAlgebra.Adjoint{T, <: AbstractVector}, b::Arr, bs...) where {T}
    return wrap(*(a, unwrap(b), bs...))
end
function *(a::LinearAlgebra.Adjoint{T, <: AbstractVector}, b::Arr, c::AbstractVector, bs...) where {T}
    return wrap(*(a, unwrap(b), unwrap(c), bs...))
end
function *(a::Number, b::Arr, bs...)
    return wrap(*(unwrap(a), unwrap(b), bs...))
end
function *(x1::Arr, x2::BasicSymbolic{VartypeT}, args...)
    return wrap(*(unwrap(x1), x2, args...))
end
function *(x1::Arr, x2::Arr, args...)
    return wrap(*(unwrap(x1), unwrap(x2), args...))
end
function *(x1::Arr, x2::AbstractMatrix, args...)
    return wrap(*(unwrap(x1), x2, args...))
end
function *(x1::Arr, x2::AbstractVector, args...)
    return wrap(*(unwrap(x1), x2, args...))
end
function *(x1::AbstractMatrix, x2::Arr, args...)
    return wrap(*(x1, unwrap(x2), args...))
end
function *(x1::Arr, x2::Arr, x3::Arr, args...)
    return wrap(*(unwrap(x1), unwrap(x2), unwrap(x3), args...))
end
function *(x1::Arr, x2::Arr, x3::Arr, x4::Arr, args...)
    return wrap(*(unwrap(x1), unwrap(x2), unwrap(x3), unwrap(x4), args...))
end

function +(x::Arr, args...)
    return +(unwrap(x), args...)
end
function +(x::Arr, y::AbstractArray, args...)
    return +(unwrap(x), y, args...)
end
function +(x1::Arr, x2::Arr, args...)
    return +(unwrap(x1), unwrap(x2), args...)
end

function +(a::SymbolicUtils.PolyadicNumericOpFirstArgT, b::Arr, bs...)
    return +(a, unwrap(b), bs...)
end

#################### MAP-REDUCE ################

SymbolicUtils.@map_methods Arr unwrap wrap
SymbolicUtils.@mapreduce_methods Arr unwrap wrap
