import Base: getindex
inner_unwrap(x) = x isa AbstractArray ? unwrap.(x) : x

##### getindex #####
@inline _indexed_ndims(idx, idxs...) = SymbolicUtils._indexed_ndims(typeof(idx)) + _indexed_ndims(idxs...)
@inline _indexed_ndims() = 0

function wrapper_fn_from_idxs(x::Arr{T, N}, idxs...) where {T, N}
    any(i -> i isa Union{Num, BasicSymbolic{VartypeT}, Arr}, idxs) && return identity
    
    nd = _indexed_ndims(idxs...)
    return nd == 0 ? is_wrapper_type(T) ? T : identity : Arr{T, nd}
end
wrapper_fn_from_idxs(x::Arr{T, N}, idx::SymbolicUtils.StableIndex{Int}) where {T, N} = T
# Wrapped array should wrap the elements too
function Base.getindex(x::Arr{T, N}, idx::CartesianIndex{N}) where {T, N}
    if is_wrapper_type(T)
        T(unwrap(x)[idx])
    else
        unwrap(x)[idx]
    end
end
function Base.getindex(x::Arr, idx...)
    wrapper_fn_from_idxs(x, idx...)(unwrap(x)[idx...])
end
const SymIdxT = Union{Num, BasicSymbolic{VartypeT}}
function Base.getindex(x::Arr, idx::SymIdxT, idxs...)
    wrapper_fn_from_idxs(x, unwrap(idx), idxs...)(unwrap(x)[idx, idxs...])
end
function Base.getindex(x::Arr, i1, idx::SymIdxT, idxs...)
    wrapper_fn_from_idxs(x, i1, unwrap(idx), idxs...)(unwrap(x)[i1, idx, idxs...])
end
function Base.getindex(x::Arr, i1::SymIdxT, idx::SymIdxT, idxs...)
    wrapper_fn_from_idxs(x, unwrap(i1), unwrap(idx), idxs...)(unwrap(x)[i1, idx, idxs...])
end
function Base.getindex(x::Arr, i1, i2, idx::SymIdxT, idxs...)
    wrapper_fn_from_idxs(x, i1, i2, unwrap(idx), idxs...)(unwrap(x)[i1, i2, idx, idxs...])
end
function Base.getindex(x::Arr, i1, i2::SymIdxT, idx::SymIdxT, idxs...)
    wrapper_fn_from_idxs(x, i1, unwrap(i2), unwrap(idx), idxs...)(unwrap(x)[i1, i2, idx, idxs...])
end
function Base.getindex(x::Arr, i1::SymIdxT, i2::SymIdxT, idx::SymIdxT, idxs...)
    wrapper_fn_from_idxs(x, unwrap(i1), unwrap(i2), unwrap(idx), idxs...)(unwrap(x)[i1, i2, idx, idxs...])
end

import Base: +, -, *

#### Broadcast ####

# On wrapper:
struct SymWrapBroadcast <: Broadcast.BroadcastStyle end

Broadcast.BroadcastStyle(::Type{T}) where {T <: Arr} = SymWrapBroadcast()

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

const PolyadicT = Union{AbstractArray{<:Number}, Number}

function *(x::Arr, args::PolyadicT...)
    return *(unwrap(x), args...)
end
function *(x::Arr, y::Number, args::PolyadicT...)
    return *(unwrap(x), unwrap(y), args...)
end
function *(a::PolyadicT, b::Arr, bs::PolyadicT...)
    return *(a, unwrap(b), bs...)
end
function *(a::LinearAlgebra.Adjoint{T, <: AbstractVector}, b::Arr, bs::PolyadicT...) where {T <: Number}
    return *(a, unwrap(b), bs...)
end
function *(a::LinearAlgebra.Adjoint{T, <: AbstractVector}, b::Arr, c::AbstractVector, bs::PolyadicT...) where {T <: Number}
    return *(a, unwrap(b), unwrap(c), bs...)
end
function *(a::Number, b::Arr, bs::PolyadicT...)
    return *(unwrap(a), unwrap(b), bs...)
end
function *(x1::Arr, x2::BasicSymbolic{VartypeT}, args::PolyadicT...)
    return *(unwrap(x1), x2, args...)
end
function *(x1::Arr, x2::Arr, args::PolyadicT...)
    return *(unwrap(x1), unwrap(x2), args...)
end
function *(x1::Arr, x2::AbstractMatrix, args::PolyadicT...)
    return *(unwrap(x1), x2, args...)
end
function *(x1::Arr, x2::AbstractVector, args::PolyadicT...)
    return *(unwrap(x1), x2, args...)
end
function *(x1::AbstractMatrix, x2::Arr, args::PolyadicT...)
    return *(x1, unwrap(x2), args...)
end
function *(x1::Arr, x2::Arr, x3::Arr, args::PolyadicT...)
    return *(unwrap(x1), unwrap(x2), unwrap(x3), args...)
end
function *(x1::Arr, x2::Arr, x3::Arr, x4::Arr, args::PolyadicT...)
    return *(unwrap(x1), unwrap(x2), unwrap(x3), unwrap(x4), args...)
end

function +(x::Arr, args::AbstractArray...)
    return +(unwrap(x), args...)
end
function +(x::Arr, y::AbstractArray, args::AbstractArray...)
    return +(unwrap(x), y, args...)
end
function +(x1::Arr, x2::Arr, args::AbstractArray...)
    return +(unwrap(x1), unwrap(x2), args...)
end

for T1 in [Arr, AbstractArray], T2 in [Arr, AbstractArray]
    T1 == T2 == AbstractArray && continue
    @eval Base.:(\)(x1::$T1{Num, 1}, x2::$T2{Num, 1}) = Num(unwrap(x1) \ unwrap(x2))
    @eval Base.:(\)(x1::$T1{Num, 1}, x2::$T2{Num, 2}) = Arr{Num, 2}(unwrap(x1) \ unwrap(x2))
    @eval Base.:(\)(x1::$T1{Num, 2}, x2::$T2{Num, 1}) = Arr{Num, 1}(unwrap(x1) \ unwrap(x2))
    @eval Base.:(\)(x1::$T1{Num, 2}, x2::$T2{Num, 2}) = Arr{Num, 2}(unwrap(x1) \ unwrap(x2))

    @eval Base.:(/)(x1::$T1{Num, 1}, x2::$T2{Num, 1}) = Arr{Num, 2}(unwrap(x1) / unwrap(x2))
    @eval Base.:(/)(x1::$T1{Num, 1}, x2::$T2{Num, 2}) = Arr{Num, 2}(unwrap(x1) / unwrap(x2))
    @eval Base.:(/)(x1::$T1{Num, 2}, x2::$T2{Num, 2}) = Arr{Num, 2}(unwrap(x1) / unwrap(x2))
end

Base.:(/)(x1::Num, x2::Arr{Num, 1}) = Arr{Num, 2}(unwrap(x1) / unwrap(x2))

#################### MAP-REDUCE ################

SymbolicUtils.@map_methods Arr unwrap wrap
SymbolicUtils.@mapreduce_methods Arr unwrap wrap

#################### REGISTRATION ################

@register_array_symbolic LinearAlgebra.triu(x::AbstractArray) begin
    size = size(x)
    eltype = eltype(x)
end

@register_array_symbolic LinearAlgebra.tril(x::AbstractArray) begin
    size = size(x)
    eltype = eltype(x)
end

@register_array_symbolic LinearAlgebra.normalize(x::AbstractArray) begin
    size = size(x)
    eltype = eltype(x)
end

@register_symbolic LinearAlgebra.tr(x::AbstractMatrix)
