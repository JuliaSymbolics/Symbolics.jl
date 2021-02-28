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
elt(s::SymArray) = elt(symtype(s))
elt(::Type{<:AbstractArray{T}}) where {T} = T
elt(::Type{<:AbstractArray}) = nothing

nd(s::SymArray) = nd(symtype(s))
nd(::Type{<:AbstractArray{<:Any, N}}) where {N} = N
nd(::Type{<:AbstractArray}) = nothing

shape(s::SymArray) = getmetadata(s, ArrayShapeCtx)

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
function promote_symtype(::typeof(getindex),
                         A::Type{<:AbstractArray},
                         idx...)
    D = count(x->x <: Number, idx)
    @maybe T=elt(A) begin
        @maybe N=nd(A) return N-D == 0 ? T : slicetype(A){T,N-D}
        return slicetype(A){T}
    end

    @maybe N=nd(A) return N-D == 0 ? T : slicetype(A){T, N-D} where T

    return slicetype(A)
end

function Base.getindex(x::SymArray, idx...)
    if any(i-> i isa Symbolic, idx)
        Term(getindex, [x, idx...])
    else
        shp = @maybe s=shape(x) s[idx...]
        setmetadata(Term(getindex, [x, idx...]),
                    ArrayShapeCtx,
                    shp)
    end
end

function Base.getindex(x::Symbolic{T}, idx::Int...) where {T<:AbstractArray}
    Term{eltype(T)}(getindex, x, [idx...])
end

# basic

# these methods are not symbolic but work if we know this info.
import Base: eltype, length, ndims, size, axes, eachindex

function eltype(A::SymArray)
    @maybe T=elt(A) return T
    error("eltype of $A not known")
end

function length(A::SymArray)
    @maybe s=shape(A) return length(s)
    error("length of $A not known")
end

function ndims(A::SymArray)
    @maybe n=nd(A) return n
    @maybe s=shape(A) return length(s.axes)
    error("ndims of $A not known")
end

function size(A::SymArray)
    @maybe s=shape(A) return length.(s.axes)
    error("size of $A not known")
end

function axes(A::SymArray)
    @maybe s=shape(A) return s.axes
    error("axes of $A not known")
end


function axes(A::SymArray, i)
    @maybe s=shape(A) begin
        @show s.axes
        return i <= length(s.axes) ? s.axes[i] : Base.OneTo(1)
    end
    error("axes of $A not known")
end

function eachindex(A::SymArray)
    @maybe s=shape(A) return CartesianIndices(s.axes)
    error("eachindex of $A not known")
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

Base.length(a::ArrayShape) = prod(map(length, axes(a)))
