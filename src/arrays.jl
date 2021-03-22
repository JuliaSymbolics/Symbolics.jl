const SymArray = Symbolic{<:AbstractArray}

struct ArrayShapeCtx end

# AbstractArray deosn't work
# indexing Matrix should give Array not AbstractArray
# length does not work

# Array interface, assumes that s.metadata is an ArrayShape, see below
# TODO: if shape is not known these should return Symbolic results

maybe(f, x) = x === nothing ? nothing : f(x)

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
function shape(s)
    ArrayShape(axes(s))
end

macro maybe(args...)
    f = args[end]
    vars = args[1:end-1]
    names = [(@assert v.head == :(=); v.args[1]) for v in vars]
    quote
        $(vars...)
        if !any(isnothing, ($(names...),))
            $f
        end
        # nothing otherwise
    end |> esc
end

slicetype(::Type{<:Array}) = Array
slicetype(::Type) = AbstractArray

struct Unknown end

propagate_atype(f, args...) = AbstractArray
propagate_eltype(f, args...) = Unknown()
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

propagate_atype(::typeof(getindex), x, idx...) = slicetype(symtype(x))
function propagate_ndims(::typeof(getindex), x, idx...)
    maybe(getndims(x)) do N
        N - count(x->symtype(x) <: Integer, idx)
    end
end

function propagate_shape(::typeof(getindex), x, idx...)
    maybe(shape(x)) do s
        s[idx...]
    end
end

function Base.getindex(x::SymArray, idx...)
    arrterm(getindex, x, idx...)
end

function Base.getindex(x::Symbolic{T}, idx::Int...) where {T<:AbstractArray}
    Term{eltype(T)}(getindex, x, [idx...])
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
    return length(s)
end

function ndims(A::SymArray)
    n = getndims(A)
    if n === Unknown()
        s = shape(A)
        if s === Unknown()
            error("ndims of $A not known")
        end
        return length(s.axes)
    end
    return n
end

function size(A::SymArray)
    s = shape(A)
    s === Unknown() && error("size of $A not known")
    return length.(s.axes)
end

function size(A::SymArray, i::Integer)
    @assert(i > 0)
    i > ndims(A) ? 1 : size(A)[i]
end

function axes(A::SymArray)
    s = shape(A)
    s === Unknown() && error("axes of $A not known")
    return s.axes
end


function axes(A::SymArray, i)
    s = shape(A)
    s === Unknown() && error("axes of $A not known")
    return i <= length(s.axes) ? s.axes[i] : Base.OneTo(1)
end

function eachindex(A::SymArray)
    s = shape(A)
    s === Unknown() && error("eachindex of $A not known")
    return CartesianIndices(s.axes)
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

# ArrayShape
# Note: implement this as if it's an array
# the idea is it needs to be usable both during construction
# and simplification

struct ArrayShape
    axes::Tuple
end

Base.axes(a::ArrayShape) = a.axes

function Base.getindex(a::ArrayShape, idx...)
    axes = a.axes
    idx1 = to_indices(CartesianIndices(axes), axes, idx)
    newaxes = ([1:length(x) for x in idx1 if !(x isa Number)]...,)
    newshape = ArrayShape(newaxes)
end

Base.ndims(a::ArrayShape) = length(a.axes)

Base.length(a::ArrayShape) = prod(map(length, axes(a)))
