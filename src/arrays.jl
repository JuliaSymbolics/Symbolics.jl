using SymbolicUtils
using StaticArrays
import Base: eltype, ndims, size

const SymArray = Symbolic{<:AbstractArray}

### Store Shape as a metadata in Term{<:AbstractArray} objects
struct ArrayShapeCtx end


#=
 There are 2 types of array terms:
   `ArrayOp{T<:AbstractArray}` and `Term{<:AbstractArray}`

 - ArrayOp represents a Einstein-notation-inspired array operation.
   it contains a field `term` which is a `Term` that represents the
   operation that resulted in the `ArrayOp`.
   I.e. will be `A*b` for the operation `(i,) => A[i,j] * b[j]` for example.
   It can be `nothing` if not known.
   - calling `shape` on an `ArrayOp` will return the shape of the array or `Unknown()`
   - do not rely on the `symtype` or `shape` information of the `.term` when looking at an `ArrayOp`.
     call `shape`, `symtype` and `ndims` directly on the `ArrayOp`.
 - `Term{<:AbstractArray}`
   - calling `shape` on it will return the shape of the array or `Unknown()`, and uses
     the `ArrayShapeCtx` metadata context to store this.

 The Array type parameter must contain the dimension.
=#

#### ArrayOp ####

"""
    ArrayOp(output_idx, expr, reduce)

A tensor expression where `output_idx` are the output indices
`expr`, is the tensor expression and `reduce` is the function
used to reduce over contracted dimensions.
"""
struct ArrayOp{T<:AbstractArray}
    output_idx # output indices
    expr       # Used in pattern matching
               # Useful to infer eltype
    reduce
    term::Ref{Any} # may or may not exist
end

symtype(a::ArrayOp{T}) where {T} = T
istree(a::ArrayOp) = true
operation(a::ArrayOp) = typeof(a)
arguments(a::ArrayOp) = [a.output_idx, a.expr, a.term, a.reduce]

macro arrayop(call, output_idx, expr, reduce=+)
    @assert output_idx.head == :tuple
    oidxs = filter(x->x isa Symbol, output_idx.args)
    iidxs = find_indices(expr)
    idxs  = union(oidxs, iidxs)
    fbody = call2term(deepcopy(expr))
    oftype(x,T) = :($x::$T)
    aop = gensym("aop")
    quote
        let
            @syms $(map(x->oftype(x, Int), idxs)...)

            $aop = $ArrayOp($output_idx,
                     $fbody,
                     $reduce,
                     $(call2term(call)))

            $aop.term[] = adjust_metadata($(call2term(call)))
            $aop
        end
    end |> esc
end

### Propagate ###
#
## Shape ##

function shape(aop::ArrayOp)
    output_idx = aop.output_idx
    expr = aop.expr

    matches = idx_to_axes(expr)
    for (sym, ms) in matches
        @assert !isempty(ms) "dimension of $sym is unknown"
        to_check = filter(m->!isnothing(shape(m.A)), ms)
        # Only check known dimensions. It may be "known symbolically"
        isempty(to_check) && continue
        reference = axes(first(to_check).A, first(to_check).dim)
        for i in 2:length(ms)
            m = ms[i]
            s=shape(m.A)
            if s !== Unknown()
                if !isequal(axes(m.A, m.dim), reference)
                    "expected axes($(m.A), $(m.dim)) = $(reference)" |> DimensionMismatch |> throw
                end
            end
        end
    end

    map(output_idx) do i
        mi = matches[i]
        @assert !isempty(mi)
        get(first(mi))
    end
end

## Eltype ##

# TODO: have fallback
function eltype(aop::ArrayOp)
    symtype(aop.expr)
end

## Ndims ##
function ndims(aop::ArrayOp)
    length(aop.output_idx)
end


### Utils ###


# turn `f(x...)` into `term(f, x...)`
#
function call2term(expr, arrs=[])
    !(expr isa Expr) && return expr
    if expr.head == :call
        return Expr(:call, term, map(call2term, expr.args)...)
    elseif expr.head == Symbol("'")
        return Expr(:call, term, adjoint, map(call2term, expr.args)...)
    end

    return Expr(expr.head, map(call2term, expr.args)...)
end

# Find all symbolic indices in expr
function find_indices(expr, idxs=[])
    !(expr isa Expr) && return idxs
    if expr.head == :ref
        return append!(idxs, filter(x->x isa Symbol, expr.args[2:end]))
    elseif expr.head == :call && expr.args[1] == :getindex || expr.args[1] == getindex
        return append!(idxs, filter(x->x isa Symbol, expr.args[3:end]))
    else
        foreach(x->find_indices(x, idxs), expr.args)
        return idxs
    end
end

struct AxisOf
    A
    dim
end

Base.get(a::AxisOf) = axes(a.A, a.dim)

function idx_to_axes(expr, dict=Dict{Sym, Vector}())
    if istree(expr)
        if operation(expr) === (getindex)
            args = arguments(expr)
            for (axis, sym) in enumerate(@views args[2:end])
                !(sym isa Sym) && continue
                axesvec = Base.get!(() -> [], dict, sym)
                push!(axesvec, AxisOf(car(args), axis))
            end
        else
            foreach(ex->idx_to_axes(ex, dict), arguments(expr))
        end
    end
    dict
end


#### Term{<:AbstractArray}
#

"""
    arrterm(f, args...; arrayop=nothing)

Create a term of `Term{<: AbstractArray}` which
is the representation of `f(args...)`.

- Calls `propagate_atype(f, args...)` to determine the
  container type, i.e. `Array` or `StaticArray` etc.
- Calls `propagate_eltype(f, args...)` to determine the
  output element type.
- Calls `propagate_ndims(f, args...)` to determine the
  output dimension.
- Calls `propagate_shape(f, args...)` to determine the
  output array shape.

`propagate_shape`, `propagate_atype`, `propagate_eltype` may
return `Unknown()` to say that the output cannot be determined

But `propagate_ndims` must work and return a non-negative integer.
"""
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

"""
    shape(s::Any)

Returns `axes(s)` or throws.
"""
shape(s) = axes(s)

"""
    shape(s::SymArray)

Extract the shape metadata from a SymArray.
If not known, returns `Unknown()`
"""
function shape(s::SymArray)
    if hasmetadata(s, ArrayShapeCtx)
        getmetadata(s, ArrayShapeCtx)
    else
        Unknown()
    end
end

## `propagate_` interface:
#  used in the `arrterm` construction.

atype(::Type{<:Array}) = Array
atype(::Type{<:SArray}) = SArray
atype(::Type) = AbstractArray

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

function propagate_ndims(f, args...)
    error("Could not determine the output dimension of $f$args")
end

function propagate_shape(f, args...)
    error("Don't know how to propagate shape for $f$args")
end

################# Base array functions
#

const ArrayLike = Union{ArrayOp, SymArray}

# basic
# these methods are not symbolic but work if we know this info.
import Base: eltype, length, ndims, size, axes, eachindex

geteltype(s::ArrayLike) = geteltype(symtype(s))
geteltype(::Type{<:AbstractArray{T}}) where {T} = T
geteltype(::Type{<:AbstractArray}) = Unknown()

ndims(s::ArrayLike) = ndims(symtype(s))
ndims(T::Type{<:AbstractArray}) = ndims(T)

function eltype(A::ArrayLike)
    T = geteltype(A)
    T === Unknown() && error("eltype of $A not known")
    return T
end

function length(A::ArrayLike)
    s = shape(A)
    s === Unknown() && error("length of $A not known")
    return prod(length, s)
end

function size(A::ArrayLike)
    s = shape(A)
    s === Unknown() && error("size of $A not known")
    return length.(s)
end

function size(A::ArrayLike, i::Integer)
    @assert(i > 0)
    i > ndims(A) ? 1 : size(A)[i]
end

function axes(A::ArrayLike)
    s = shape(A)
    s === Unknown() && error("axes of $A not known")
    return s
end


function axes(A::ArrayLike, i)
    s = shape(A)
    s === Unknown() && error("axes of $A not known")
    return i <= length(s) ? s[i] : Base.OneTo(1)
end

function eachindex(A::ArrayLike)
    s = shape(A)
    s === Unknown() && error("eachindex of $A not known")
    return CartesianIndices(s)
end
