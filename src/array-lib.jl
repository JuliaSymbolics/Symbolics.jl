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
    wrapper_fn_from_idxs(x, idx...)(unwrap(x)[unwrap.(idx)...])
end
const SymIdxT = Union{Num, BasicSymbolic{VartypeT}}
function Base.getindex(x::Arr, idx::SymIdxT, idxs...)
    wrapper_fn_from_idxs(x, unwrap(idx), idxs...)(unwrap(x)[unwrap(idx), unwrap.(idxs)...])
end
function Base.getindex(x::Arr, i1, idx::SymIdxT, idxs...)
    wrapper_fn_from_idxs(x, i1, unwrap(idx), idxs...)(unwrap(x)[(i1), unwrap(idx), unwrap.(idxs)...])
end
function Base.getindex(x::Arr, i1::SymIdxT, idx::SymIdxT, idxs...)
    wrapper_fn_from_idxs(x, unwrap(i1), unwrap(idx), idxs...)(unwrap(x)[unwrap(i1), unwrap(idx), unwrap.(idxs)...])
end
function Base.getindex(x::Arr, i1, i2, idx::SymIdxT, idxs...)
    wrapper_fn_from_idxs(x, i1, i2, unwrap(idx), idxs...)(unwrap(x)[unwrap(i1), unwrap(i2), unwrap(idx), unwrap.(idxs)...])
end
function Base.getindex(x::Arr, i1::SymIdxT, i2, idx::SymIdxT, idxs...)
    return wrapper_fn_from_idxs(x, unwrap(i1), i2, unwrap(idx), idxs...)(unwrap(x)[unwrap(i1), unwrap(i2), unwrap(idx), unwrap.(idxs)...])
end
function Base.getindex(x::Arr, i1, i2::SymIdxT, idx::SymIdxT, idxs...)
    wrapper_fn_from_idxs(x, i1, unwrap(i2), unwrap(idx), idxs...)(unwrap(x)[unwrap(i1), unwrap(i2), unwrap(idx), unwrap.(idxs)...])
end
function Base.getindex(x::Arr, i1::SymIdxT, i2::SymIdxT, idx::SymIdxT, idxs...)
    wrapper_fn_from_idxs(x, unwrap(i1), unwrap(i2), unwrap(idx), idxs...)(unwrap(x)[unwrap(i1), unwrap(i2), unwrap(idx), unwrap.(idxs)...])
end

import Base: +, -, *

#### Broadcast ####

# On wrapper:
struct SymWrapBroadcast <: Broadcast.BroadcastStyle end

Broadcast.BroadcastStyle(::Type{T}) where {T <: Arr} = SymWrapBroadcast()

Broadcast.BroadcastStyle(::SymWrapBroadcast,
    ::Broadcast.BroadcastStyle) = SymWrapBroadcast()
Broadcast.BroadcastStyle(::SymWrapBroadcast, ::Broadcast.Unknown) = SymWrapBroadcast()
Broadcast.BroadcastStyle(::SymbolicUtils.SymBroadcast,
    ::SymWrapBroadcast) = Broadcast.Unknown()

unwrap_broadcasts(head, args...) = (unwrap_broadcast(head), unwrap_broadcasts(args...)...)
unwrap_broadcasts() = ()
unwrap_broadcast(x) = unwrap(x)
function unwrap_broadcast(bc::Broadcast.Broadcasted{SymWrapBroadcast})
    Broadcast.Broadcasted{SymbolicUtils.SymBroadcast{VartypeT}}(bc.f, unwrap_broadcasts(bc.args...), bc.axes)
end

function Base.copy(bc::Broadcast.Broadcasted{SymWrapBroadcast})
    return wrap(copy(unwrap_broadcast(bc)))
end

#################### POLYADIC ################

*(x::Arr, y::Number) = *(unwrap(x), unwrap(y))
*(x::Number, y::Arr) = *(unwrap(x), unwrap(y))
*(x::Arr, y::BasicSymbolic{VartypeT}) = *(unwrap(x), y)
*(x::Arr, y::Arr) = *(unwrap(x), unwrap(y))
*(x::Arr, y::AbstractMatrix) = *(unwrap(x), y)
*(x::Arr, y::AbstractVector) = *(unwrap(x), y)
*(x::AbstractMatrix, y::Arr) = *(x, unwrap(y))
*(x::LinearAlgebra.Diagonal, y::Arr) = *(x, unwrap(y))
for T in (LinearAlgebra.Adjoint, LinearAlgebra.Transpose), N in (1, 2)
    @eval *(
        x::$T{<:Any, <:AbstractVector}, y::Arr{<:Any, $N}
    ) = *(x, unwrap(y))
end
function *(
        x::Union{
            LinearAlgebra.Adjoint{<:Any, <:AbstractVector},
            LinearAlgebra.Transpose{<:Any, <:AbstractVector},
        }, y::Arr{<:Any, 1}
    )
    return *(x, unwrap(y))
end
function *(x::LinearAlgebra.Adjoint{T, <:AbstractVector}, y::Arr{S, 1}) where {
        T <: Number, S <: Number,
    }
    return *(x, unwrap(y))
end
function *(x::LinearAlgebra.Transpose{T, <:AbstractVector}, y::Arr{T, 1}) where {T <: Real}
    return *(x, unwrap(y))
end
# StaticArrays dispatches on this non-public abstract storage type, so exact
# intersections with its methods cannot be expressed through the public aliases.
const _StaticArray = StaticArraysCore.StaticArray
function *(
        x::_StaticArray{Tuple{N, M}, T, 2}, y::Arr{S, 1}
    ) where {N, M, T, S}
    return *(x, unwrap(y))
end
function *(
        x::Arr, y::Union{
            LinearAlgebra.Adjoint{<:Any, <:AbstractVector},
            LinearAlgebra.Transpose{<:Any, <:AbstractVector},
        }
    )
    return *(unwrap(x), y)
end
function *(
        x::LinearAlgebra.Adjoint{<:Number, <:AbstractVector}, y::Arr,
        z::AbstractVector{<:Number}
    )
    return *(x, unwrap(y), unwrap(z))
end
function *(
        x::LinearAlgebra.Transpose{<:Real, <:AbstractVector}, y::Arr,
        z::AbstractVector{<:Real}
    )
    return *(x, unwrap(y), unwrap(z))
end
function *(
        x::Union{Real, Complex}, y::Arr{T, 2}, z::Arr{S, 1}
    ) where {T <: Union{Real, Complex}, S <: Union{Real, Complex}}
    return *(unwrap(x), unwrap(y), unwrap(z))
end
*(x::Arr, y::Arr, z::Number) = *(unwrap(x), unwrap(y), unwrap(z))

function +(x::Arr, args::AbstractArray...)
    return +(unwrap(x), args...)
end
function +(x::Arr, y::AbstractArray, args::AbstractArray...)
    return +(unwrap(x), y, args...)
end
function +(x1::Arr, x2::Arr, args::AbstractArray...)
    return +(unwrap(x1), unwrap(x2), args...)
end
function +(x::Arr, y::_StaticArray)
    return +(unwrap(x), y)
end

for T1 in [Arr, AbstractArray], T2 in [Arr, AbstractArray]
    T1 == T2 == AbstractArray && continue
    @eval Base.:(\)(x1::$T1{Num, 1}, x2::$T2{Num, 1}) = Num(unwrap(x1) \ unwrap(x2))
    @eval Base.:(\)(x1::$T1{Num, 1}, x2::$T2{Num, 2}) = Arr{Num, 2}(unwrap(x1) \ unwrap(x2))

    @eval Base.:(/)(x1::$T1{Num, 1}, x2::$T2{Num, 1}) = Arr{Num, 2}(unwrap(x1) / unwrap(x2))
    @eval Base.:(/)(x1::$T1{Num, 1}, x2::$T2{Num, 2}) = Arr{Num, 2}(unwrap(x1) / unwrap(x2))
    @eval Base.:(/)(x1::$T1{Num, 2}, x2::$T2{Num, 2}) = Arr{Num, 2}(unwrap(x1) / unwrap(x2))
end

Base.:(/)(x1::Num, x2::Arr{Num, 1}) = Arr{Num, 2}(unwrap(x1) / unwrap(x2))

const _TriangularNum = Union{
    LinearAlgebra.LowerTriangular{Num}, LinearAlgebra.UpperTriangular{Num},
}
const _UnitTriangularNum = Union{
    LinearAlgebra.UnitLowerTriangular{Num}, LinearAlgebra.UnitUpperTriangular{Num},
}
for T in (LinearAlgebra.Diagonal{Num}, _TriangularNum, _UnitTriangularNum), N in (1, 2)
    @eval Base.:(/)(A::Arr{Num, $N}, B::$T) = Arr{Num, 2}(unwrap(A) / unwrap(B))
end
for T in (
        LinearAlgebra.Bidiagonal{Num},
        LinearAlgebra.Adjoint{Num, <:LinearAlgebra.Bidiagonal},
        LinearAlgebra.Transpose{Num, <:LinearAlgebra.Bidiagonal},
    )
    @eval Base.:(/)(A::Arr{Num, 2}, B::$T) = Arr{Num, 2}(unwrap(A) / unwrap(B))
end
for T in (LinearAlgebra.Adjoint, LinearAlgebra.Transpose)
    @eval function Base.:(/)(
            A::$T{Num, <:AbstractVector}, B::Arr{Num, 2}
        )
        return Arr{Num, 2}(unwrap(A) / unwrap(B))
    end
end

Base.exp(m::Matrix{Num}) = Arr{Num, 2}(exp(SConst(m)))
Base.exp(m::Matrix{Complex{Num}}) = Arr{Complex{Num}, 2}(exp(SConst(m)))

function LinearAlgebra.transpose(x::Arr{T, N}) where {T, N}
    return Arr{T, 2}(LinearAlgebra.transpose(unwrap(x)))
end

#################### MAP-REDUCE ################

for Tf in (BasicSymbolic, Any)
    @eval begin
        function Base.map(f::$Tf, x::Arr)
            return wrap(map(f, unwrap(x)))
        end
        function Base.map(f::$Tf, x::Arr, xs::Arr...)
            return wrap(map(f, unwrap(x), unwrap.(xs)...))
        end
        function Base.map(f::$Tf, x::AbstractArray, y::Arr, ys::Arr...)
            return wrap(map(f, x, unwrap(y), unwrap.(ys)...))
        end
        function Base.map(f::$Tf, x::_StaticArray, y::Arr, ys::Arr...)
            return wrap(map(f, x, unwrap(y), unwrap.(ys)...))
        end
        function Base.map(f::$Tf, x::BasicSymbolic, y::Arr, ys::Arr...)
            return wrap(map(f, x, unwrap(y), unwrap.(ys)...))
        end
    end
end
for Tf in (BasicSymbolic, Any), Tr in (BasicSymbolic, Any)
    @eval begin
        function Base.mapreduce(f::$Tf, op::$Tr, x::Arr; kw...)
            return wrap(mapreduce(f, op, unwrap(x); kw...))
        end
        function Base.mapreduce(f::$Tf, op::$Tr, x::Arr, xs::Arr...; kw...)
            return wrap(mapreduce(f, op, unwrap(x), unwrap.(xs)...; kw...))
        end
        function Base.mapreduce(
                f::$Tf, op::$Tr, x::AbstractArray, y::Arr, ys::Arr...; kw...
            )
            return wrap(mapreduce(f, op, x, unwrap(y), unwrap.(ys)...; kw...))
        end
        function Base.mapreduce(
                f::$Tf, op::$Tr, x::BasicSymbolic, y::Arr, ys::Arr...; kw...
            )
            return wrap(mapreduce(f, op, x, unwrap(y), unwrap.(ys)...; kw...))
        end
    end
end

function LinearAlgebra.dot(x::Arr{T}, y::Arr{T}) where {T}
    T(LinearAlgebra.dot(unwrap(x), unwrap(y)))
end

#################### REGISTRATION ################

@register_array_symbolic LinearAlgebra.triu(x::AbstractArray) begin
    size = size(x)
    eltype = eltype(x)
    ndims = ndims(x)
end

@register_array_symbolic LinearAlgebra.tril(x::AbstractArray) begin
    size = size(x)
    eltype = eltype(x)
    ndims = ndims(x)
end

@register_array_symbolic LinearAlgebra.normalize(x::AbstractArray) begin
    size = size(x)
    eltype = eltype(x)
    ndims = ndims(x)
end

@register_symbolic LinearAlgebra.tr(x::AbstractMatrix)

# Eigenvalues of a symbolic matrix have no closed form, so `eigmax`/`eigmin`
# build unevaluated symbolic terms (like `tr`) rather than computing numerically;
# without this they fall through to `eigvals!`, which errors on a symbolic matrix.
@register_symbolic LinearAlgebra.eigmax(x::AbstractMatrix)
@register_symbolic LinearAlgebra.eigmin(x::AbstractMatrix)
