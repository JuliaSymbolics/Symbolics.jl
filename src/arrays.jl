const SymArray = Symbolic{<:AbstractArray}
using StaticArrays

struct ArrayShapeCtx end

# AbstractArray deosn't work
# indexing Matrix should give Array not AbstractArray
# length does not work

# Array interface, assumes that s.metadata is an ArrayShape, see below
# TODO: if shape is not known these should return Symbolic results

# allows promote_symtype to be called with the correct types
symtype(x::Union{Colon, AbstractRange}) = typeof(x)


# Partial information
geteltype(s::SymArray) = geteltype(symtype(s))
geteltype(::Type{<:AbstractArray{T}}) where {T} = T
geteltype(::Type{<:AbstractArray}) = Unknown()

getndims(s::SymArray) = getndims(symtype(s))
getndims(x) = ndims(x)
getndims(::Type{<:AbstractArray{<:Any, N}}) where {N} = N
getndims(::Type{<:AbstractArray}) = Unknown()

function shape(s::SymArray)
    if hasmetadata(s, ArrayShapeCtx)
        getmetadata(s, ArrayShapeCtx)
    else
        Unknown()
    end
end
shape(s) = axes(s)

struct Unknown end

_propagate_atype(::Type{T}, ::Type{T}) where {T} = T
_propagate_atype(::Type{<:Array}, ::Type{<:SArray}) = Array
_propagate_atype(::Type{<:SArray}, ::Type{<:Array}) = Array
_propagate_atype(::Any, ::Any) = AbstractArray
_propagate_atype(T) = T
_propagate_atype() = AbstractArray

function propagate_atype(f, args...)
    As = [atype(symtype(T))
          for T in Iterators.filter(x->x <: Symbolic{<:AbstractArray}, typeof.(args))]
    if length(As) <= 1
        _propagate_atype(As...)
    else
        foldl(_propagate_atype, As)
    end
end

function propagate_eltype(f, args...)
    As = [eltype(symtype(T))
          for T in Iterators.filter(x->symtype(x) <: AbstractArray, args)]
    promote_type(As...)
end

propagate_ndims(f, args...) = Unknown()

function propagate_shape(f, args...)
    error("Don't know how to propagate shape for $f$args")
end

function arrterm(f, args...)
    atype = propagate_atype(f, args...)
    etype = propagate_eltype(f, args...)
    nd    = propagate_ndims(f, args...)

    S = if etype === Unknown() && nd === Unknown()
        atype
    elseif etype === Unknown()
        atype{T, nd} where T
    elseif nd === Unknown()
        atype{etype, N} where N
    else
        atype{etype, nd}
    end

    setmetadata(Term{S}(f, args),
                ArrayShapeCtx,
                propagate_shape(f, args...))
end

maybe(f, x::Unknown) = Unknown()
maybe(f, x) = f(x)

function maybefoldl(f, g, xs, acc)
    for x in xs
        y = f(x)
        y === Unknown() && return Unknown()
        acc = g(acc, y)
    end
    return acc
end
atype(::Type{<:Array}) = Array
atype(::Type{<:SArray}) = SArray
atype(::Type) = AbstractArray

function propagate_ndims(::typeof(getindex), x, idx...)
    maybe(getndims(x)) do N
        N - count(x->symtype(x) <: Integer, idx)
    end
end

function propagate_shape(::typeof(getindex), x, idx...)
    axes = shape(x)
    axes === Unknown() && return Unknown()

    idx1 = to_indices(CartesianIndices(axes), axes, idx)
    ([1:length(x) for x in idx1 if !(symtype(x) <: Number)]...,)
end

function Base.getindex(x::SymArray, idx...)
    if all(i->symtype(i) <: Integer, idx)
        Term{eltype(symtype(x))}(getindex, [x, idx...])
    else
        arrterm(getindex, x, idx...)
    end
end

# basic
# these methods are not symbolic but work if we know this info.
import Base: eltype, length, ndims, size, axes, eachindex

function eltype(A::SymArray)
    T = geteltype(A)
    T === Unknown() && error("eltype of $A not known")
    return T
end

function length(A::SymArray)
    s = shape(A)
    s === Unknown() && error("length of $A not known")
    return prod(length, s)
end

function ndims(A::SymArray)
    n = getndims(A)
    if n === Unknown()
        s = shape(A)
        if s === Unknown()
            error("ndims of $A not known")
        end
        return length(s)
    end
    return n
end

function size(A::SymArray)
    s = shape(A)
    s === Unknown() && error("size of $A not known")
    return length.(s)
end

function size(A::SymArray, i::Integer)
    @assert(i > 0)
    i > ndims(A) ? 1 : size(A)[i]
end

function axes(A::SymArray)
    s = shape(A)
    s === Unknown() && error("axes of $A not known")
    return s
end


function axes(A::SymArray, i)
    s = shape(A)
    s === Unknown() && error("axes of $A not known")
    return i <= length(s) ? s[i] : Base.OneTo(1)
end

function eachindex(A::SymArray)
    s = shape(A)
    s === Unknown() && error("eachindex of $A not known")
    return CartesianIndices(s)
end

# todo: stride?

# map, reduce
# reshape
# copy
# deepcopy
# similar
# reinterpret -- ??
# cat

# Unary arithmetic – -, +
# Binary arithmetic – -, +, *, /, \, ^
# Comparison – ==, !=, ≈ (isapprox), ≉
# Broadcast
#
