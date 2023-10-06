using SymbolicUtils
using SymbolicUtils: @capture
using StaticArrays
import Base: eltype, length, ndims, size, axes, eachindex

export @arrayop, ArrayMaker, @makearray, @setview, @setview!

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
struct ArrayOp{T<:AbstractArray} <: Symbolic{T}
    output_idx # output indices
    expr       # Used in pattern matching
               # Useful to infer eltype
    reduce
    term
    shape
    ranges::Dict{BasicSymbolic, AbstractRange} # index range each index symbol can take,
                                     # optional for each symbol
    metadata
end

function ArrayOp(T::Type, output_idx, expr, reduce, term, ranges=Dict(); metadata=nothing)
    sh = make_shape(output_idx, unwrap(expr), ranges)
    ArrayOp{T}(output_idx, unwrap(expr), reduce, term, sh, ranges, metadata)
end

function ArrayOp(a::AbstractArray)
    i = makesubscripts(ndims(a))
    # TODO: formalize symtype(::Type) then!
    ArrayOp(symtype(a), (i...,), a[i...], +, a)
end

ConstructionBase.constructorof(s::Type{<:ArrayOp{T}}) where {T} = ArrayOp{T}

shape(aop::ArrayOp) = aop.shape

const show_arrayop = Ref{Bool}(false)
function Base.show(io::IO, aop::ArrayOp)
    if istree(aop.term) && !show_arrayop[]
        show(io, aop.term)
    else
        print(io, "@arrayop")
        print(io, "(_[$(join(string.(aop.output_idx), ","))] := $(aop.expr))")
        if aop.reduce != +
            print(io, " ($(aop.reduce))")
        end

        if !isempty(aop.ranges)
            print(io, " ", join(["$k in $v" for (k, v) in aop.ranges], ", "))
        end
    end
end

Base.summary(io::IO, aop::ArrayOp) = Base.array_summary(io, aop, shape(aop))
function Base.showarg(io::IO, aop::ArrayOp, toplevel)
    show(io, aop)
    toplevel && print(io, "::", typeof(aop))
    return nothing
end

symtype(a::ArrayOp{T}) where {T} = T
istree(a::ArrayOp) = true
function operation(a::ArrayOp)
    isnothing(a.term) ? typeof(a) : operation(a.term)
end
function arguments(a::ArrayOp)
    isnothing(a.term) ? [a.output_idx, a.expr, a.reduce,
                         a.term, a.shape, a.ranges, metadata(a)] :
    arguments(a.term)
end

function similarterm(a::ArrayOp, f, args, type; metadata=nothing)
    res = f(args...)
    if res isa Symbolic && metadata !== nothing
        res = SymbolicUtils.metadata(res, metadata)
    end
    res
end

function Base.isequal(a::ArrayOp, b::ArrayOp)
    a === b && return true
    isequal(a.shape, b.shape) &&
    isequal(a.ranges, b.ranges) &&
    isequal(a.output_idx, b.output_idx) &&
    isequal(a.reduce, b.reduce) &&
    isequal(operation(a), operation(b)) &&
    isequal(a.expr, b.expr)
end

function Base.hash(a::ArrayOp, u::UInt)
    hash(a.shape, hash(a.ranges, hash(a.expr, hash(a.output_idx, hash(operation(a), u)))))
end

macro arrayop(output_idx, expr, options...)
    rs = []
    reduce = +
    call = nothing

    extra = []
    for o in options
        if isexpr(o, :call) && o.args[1] == :in
            push!(rs, :($(o.args[2]) => $(o.args[3])))
        elseif isexpr(o, :(=)) && o.args[1] == :reduce
            reduce = o.args[2]
        elseif isexpr(o, :(=)) && o.args[1] == :term
            call = o.args[2]
        else
            push!(extra, o)
        end
    end
    if length(extra) == 1
        @warn("@arrayop <call> <idx> <expr> is deprecated, use @arrayop <idx> <expr> term=<call> instead")
        call = output_idx
        output_idx = expr
        expr = extra[1]
    end
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

            expr = $fbody
            #TODO: proper Atype
            $ArrayOp(Array{$symtype(expr),
                           $(length(output_idx.args))},
                     $output_idx,
                     expr,
                     $reduce,
                     $(call2term(call)),
                     Dict($(rs...)))

        end
    end |> esc
end

const SymArray = Union{ArrayOp, Symbolic{<:AbstractArray}}
const SymMat = Union{ArrayOp{<:AbstractMatrix}, Symbolic{<:AbstractMatrix}}
const SymVec = Union{ArrayOp{<:AbstractVector}, Symbolic{<:AbstractVector}}

### Propagate ###
#
## Shape ##


function axis_in(a, b)
    first(a) >= first(b) && last(a) <= last(b)
end

function make_shape(output_idx, expr, ranges=Dict())
    matches = idx_to_axes(expr)
    for (sym, ms) in matches
        to_check = filter(m->!(shape(m.A) isa Unknown), ms)
        # Only check known dimensions. It may be "known symbolically"
        isempty(to_check) && continue
        restricted = false
        if haskey(ranges, sym)
            ref_axis = ranges[sym]
            restricted = true
        else
            ref_axis = axes(first(to_check).A, first(to_check).dim)
        end
        reference = ref_axis
        for i in (restricted ? 1 : 2):length(ms)
            m = ms[i]
            s=shape(m.A)
            if s !== Unknown()
                if restricted
                    if !axis_in(ref_axis, axes(m.A, m.dim))
                        throw(DimensionMismatch("expected $(ref_axis) to be within axes($(m.A), $(m.dim))"))
                    end
                elseif !isequal(axes(m.A, m.dim), reference)
                    throw(DimensionMismatch("expected axes($(m.A), $(m.dim)) = $(reference)"))
                end
            end
        end
    end

    sz = map(output_idx) do i
        if issym(i)
            if haskey(ranges, i)
                return axes(ranges[i], 1)
            end
            if !haskey(matches, i)
                error("There was an error processing arrayop expression $expr.\n" *
                      "Dimension of output index $i in $output_idx could not be inferred")
            end
            mi = matches[i]
            @assert !isempty(mi)
            ext = get_extents(mi)
            ext isa Unknown && return Unknown()
            return Base.OneTo(length(ext))
        elseif i isa Integer
            return Base.OneTo(1)
        end
    end
    # TODO: maybe we can remove this restriction?
    if any(x->x isa Unknown, sz)
        Unknown()
    else
        sz
    end
end


function ranges(a::ArrayOp)
    rs = Dict{BasicSymbolic, Any}()
    ax = idx_to_axes(a.expr)
    for i in keys(ax)
        if haskey(a.ranges, i)
            rs[i] = a.ranges[i]
        else
            rs[i] = ax[i] #get_extents(ax[i])
        end
    end
    return rs
end

## Eltype ##

eltype(aop::ArrayOp) = symtype(aop.expr)

## Ndims ##
function ndims(aop::ArrayOp)
    length(aop.output_idx)
end


### Utils ###


# turn `f(x...)` into `term(f, x...)`
#
function call2term(expr, arrs=[])
    !(expr isa Expr) && return :($unwrap($expr))
    if expr.head == :call
        if expr.args[1] == :(:)
            return expr
        end
        return Expr(:call, term, map(call2term, expr.args)...)
    elseif expr.head == :ref
        return Expr(:ref, call2term(expr.args[1]), expr.args[2:end]...)
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
    boundary
end

function Base.get(a::AxisOf)
    @oops shape(a.A)
    axes(a.A, a.dim)
end

function get_extents(xs)
    boundaries = map(x->x.boundary, xs)
    if all(iszero∘wrap, boundaries)
        get(first(xs))
    else
        ii = findfirst(x->issym(x) || istree(x), boundaries)
        if !isnothing(ii)
            error("Could not find the boundary from symbolic index $(xs[ii]). Please manually specify the range of indices.")
        end
        extent = get(first(xs))
        start_offset = -reduce(min, filter(x->x<0, boundaries), init=0)
        end_offset = reduce(max, filter(x->x>0, boundaries), init=0)

        (first(extent) + start_offset):(last(extent) - end_offset)
    end
end

get_extents(x::AbstractRange) = x

## Walk expr looking for symbols used in getindex expressions
# Returns a dictionary of Sym to a vector of AxisOf objects.
# The vector has as many elements as the number of times the symbol
# appears in the expr. AxisOf has three fields:
# A: the array whose indexing it appears in
# dim: The dimension of the array indexed
# boundary: how much padding is this indexing requiring, for example
#   boundary is 2 for x[i + 2], and boundary = -2 for x[i - 2]
function idx_to_axes(expr, dict=Dict{Any, Vector}(), ranges=Dict())
    if istree(expr)
        if operation(expr) === (getindex)
            args = arguments(expr)
            for (axis, idx_expr) in enumerate(@views args[2:end])
                if issym(idx_expr) || istree(idx_expr)
                    vs = get_variables(idx_expr)
                    isempty(vs) && continue
                    sym = only(get_variables(idx_expr))
                    axesvec = Base.get!(() -> [], dict, sym)
                    push!(axesvec, AxisOf(first(args), axis, idx_expr - sym))
                end
            end
        else
            idx_to_axes(operation(expr), dict)
            foreach(ex->idx_to_axes(ex, dict), arguments(expr))
        end
    end
    dict
end


#### Term{<:AbstractArray}
#

"""
    array_term(f, args...;
        container_type = propagate_atype(f, args...),
        eltype = propagate_eltype(f, args...),
        size = map(length, propagate_shape(f, args...)),
        ndims = propagate_ndims(f, args...))

Create a term of `Term{<: AbstractArray}` which
is the representation of `f(args...)`.

Default arguments:
- `container_type=propagate_atype(f, args...)` - the container type,
    i.e. `Array` or `StaticArray` etc.
- `eltype=propagate_eltype(f, args...)` - the output element type.
- `size=map(length, propagate_shape(f, args...))` -  the
  output array size. `propagate_shape` returns a tuple of index ranges.
- `ndims=propagate_ndims(f, args...)` the output dimension.

`propagate_shape`, `propagate_atype`, `propagate_eltype` may
return `Unknown()` to say that the output cannot be determined
"""
function array_term(f, args...;
        container_type = propagate_atype(f, args...),
        eltype = propagate_eltype(f, args...),
        size = Unknown(),
        ndims = size !== Unknown() ? length(size) : propagate_ndims(f, args...),
        shape = size !== Unknown() ? Tuple(map(x->1:x, size)) : propagate_shape(f, args...))

    if container_type == Unknown()
        # There's always a fallback for this
        container_type = propagate_atype(f, args...)
    end

    if eltype == Unknown()
        eltype = Base.propagate_eltype(container_type)
    end

    if ndims == Unknown()
        ndims = if shape == Unknown()
            Any
        else
            length(shape)
        end
    end
    S = container_type{eltype, ndims}
    setmetadata(Term{S}(f, Any[args...]), ArrayShapeCtx, shape)
end

"""
    shape(s::Any)

Returns `axes(s)` or Unknown().
"""
shape(s) = axes(s)

"""
    shape(s::SymArray)

Extract the shape metadata from a SymArray.
If not known, returns `Unknown()`
"""
function shape(s::Symbolic{<:AbstractArray})
    if hasmetadata(s, ArrayShapeCtx)
        getmetadata(s, ArrayShapeCtx)
    else
        Unknown()
    end
end

## `propagate_` interface:
#  used in the `array_term` construction.

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
    if propagate_shape(f, args...) == Unknown()
        error("Could not determine the output dimension of $f$args")
    else
        length(propagate_shape(f, args...))
    end
end

function propagate_shape(f, args...)
    error("Don't know how to propagate shape for $f$args")
end

### Wrapper type for dispatch

@symbolic_wrap struct Arr{T,N} <: AbstractArray{T, N}
    value
end

Base.hash(x::Arr, u::UInt) = hash(unwrap(x), u)
Base.isequal(a::Arr, b::Arr) = isequal(unwrap(a), unwrap(b))

ArrayOp(x::Arr) = unwrap(x)

function Arr(x)
    A = symtype(x)
    @assert A <: AbstractArray
    Arr{maybewrap(eltype(A)), ndims(A)}(x)
end

const ArrayLike{T,N} = Union{
    ArrayOp{AbstractArray{T,N}},
    Symbolic{AbstractArray{T,N}},
    Arr{T,N},
    SymbolicUtils.Term{AbstractArray{T, N}}
} # Like SymArray but includes Arr and Term{Arr}

unwrap(x::Arr) = x.value

maybewrap(T) = has_symwrapper(T) ? wrapper_type(T) : T
# These methods allow @wrapped methods to be more specific and not overwrite
# each other when defined both for matrix and vector
wrapper_type(::Type{<:AbstractMatrix}) = Arr{<:Any, 2}
wrapper_type(::Type{<:AbstractMatrix{T}}) where {T} = Arr{maybewrap(T), 2}

wrapper_type(::Type{<:AbstractVector}) = Arr{<:Any, 1}
wrapper_type(::Type{<:AbstractVector{T}}) where {T} = Arr{maybewrap(T), 1}

function Base.show(io::IO, arr::Arr)
    x = unwrap(arr)
    istree(x) && print(io, "(")
    print(io, unwrap(arr))
    istree(x) && print(io, ")")
    if !(shape(x) isa Unknown)
        print(io, "[", join(string.(axes(arr)), ","), "]")
    end
end
Base.show(io::IO, ::MIME"text/plain", arr::Arr) = show(io, arr)

################# Base array functions
#

# basic
# these methods are not symbolic but work if we know this info.

geteltype(s::SymArray) = geteltype(symtype(s))
geteltype(::Type{<:AbstractArray{T}}) where {T} = T
geteltype(::Type{<:AbstractArray}) = Unknown()

ndims(s::SymArray) = ndims(symtype(s))
ndims(::Type{<:Arr{<:Any, N}}) where N = N

function eltype(A::Union{Arr, SymArray})
    T = geteltype(unwrap(A))
    T === Unknown() && error("eltype of $A not known")
    return T
end

function length(A::Union{Arr, SymArray})
    s = shape(unwrap(A))
    s === Unknown() && error("length of $A not known")
    return prod(length, s)
end

function size(A::Union{Arr, SymArray})
    s = shape(unwrap(A))
    s === Unknown() && error("size of $A not known")
    return length.(s)
end

function size(A::SymArray, i::Integer)
    @assert(i > 0)
    i > ndims(A) ? 1 : size(A)[i]
end

function axes(A::Union{Arr, SymArray})
    s = shape(unwrap(A))
    s === Unknown() && error("axes of $A not known")
    return s
end


function axes(A::SymArray, i)
    s = shape(A)
    s === Unknown() && error("axes of $A not known")
    return i <= length(s) ? s[i] : Base.OneTo(1)
end

function eachindex(A::Union{Arr, SymArray})
    s = shape(unwrap(A))
    s === Unknown() && error("eachindex of $A not known")
    return CartesianIndices(s)
end

function get_variables!(vars, e::Arr, varlist=nothing)
    foreach(x -> get_variables!(vars, x, varlist), collect(e))
    vars
end


### Scalarize

scalarize(a::Array) = map(scalarize, a)
scalarize(term::Symbolic{<:AbstractArray}, idx) = term[idx...]
val2num(::Val{n}) where n = n

function replace_by_scalarizing(ex, dict)
    rule = @rule(getindex(~x, ~~i) =>
                 scalarize(~x, (map(j->substitute(j, dict), ~~i)...,)))

    simterm = (x, f, args; kws...) -> begin
        if metadata(x) !== nothing
            similarterm(x, f, args; metadata=metadata(x))
        else
            f(args...)
        end
    end

    function rewrite_operation(x)
        if istree(x) && istree(operation(x))
            f = operation(x)
            ff = replace_by_scalarizing(f, dict)
            if metadata(x) !== nothing
                similarterm(x, ff, arguments(x); metadata=metadata(x))
            else
                ff(arguments(x)...)
            end
        else
            nothing
        end
    end

    prewalk_if(x->!(x isa ArrayOp || x isa ArrayMaker),
               Rewriters.PassThrough(Chain([rewrite_operation, rule])),
              ex, simterm)
end

function prewalk_if(cond, f, t, similarterm)
    t′ = cond(t) ? f(t) : return t
    if istree(t′)
        return similarterm(t′, operation(t′),
                           map(x->prewalk_if(cond, f, x, similarterm), arguments(t′)))
    else
        return t′
    end
end


function scalarize(arr::AbstractArray, idx)
    arr[idx...]
end

function scalarize(arr, idx)
    if istree(arr)
        scalarize_op(operation(arr), arr, idx)
    else
        error("scalarize is not defined for $arr at idx=$idx")
    end
end

scalarize_op(f, arr) = arr

struct ScalarizeCache end

function scalarize_op(f, arr, idx)
    if hasmetadata(arr, ScalarizeCache) && getmetadata(arr, ScalarizeCache)[] !== nothing
        getmetadata(arr, ScalarizeCache)[][idx...]
    else
        # wrap and unwrap to call generic methods
        thing = unwrap(f(scalarize.(map(wrap, arguments(arr)))...))
        if metadata(arr) != nothing
            # forward any metadata
            try
                thing = metadata(thing, metadata(arr))
            catch err
                @warn "could not attach metadata of subexpression $arr to the scalarized form at idx"
            end
        end
        if !hasmetadata(arr, ScalarizeCache)
            arr = setmetadata(arr, ScalarizeCache, Ref{Any}(nothing))
        end
        getmetadata(arr, ScalarizeCache)[] = thing
        thing[idx...]
    end
end

@wrapped function Base.:(\)(A::AbstractMatrix, b::AbstractVecOrMat)
    t = array_term(\, A, b)
    setmetadata(t, ScalarizeCache, Ref{Any}(nothing))
end

@wrapped function Base.inv(A::AbstractMatrix)
    t = array_term(inv, A)
    setmetadata(t, ScalarizeCache, Ref{Any}(nothing))
end

_det(x, lp) = det(x, laplace=lp)

function scalarize_op(f::typeof(_det), arr)
    unwrap(det(map(wrap, collect(arguments(arr)[1])), laplace=arguments(arr)[2]))
end

@wrapped function LinearAlgebra.det(x::AbstractMatrix; laplace=true)
    Term{eltype(x)}(_det, [x, laplace])
end


# A * x = b
# A ∈ R^(m x n) x ∈ R^(n, k) = b ∈ R^(m, k)
propagate_ndims(::typeof(\), A, b) = ndims(b)
propagate_ndims(::typeof(inv), A) = ndims(A)

# A(m,k) * B(k,n) = C(m,n)
# A(m,k) \ C(m,n)  = B(k,n)
function propagate_shape(::typeof(\), A, b)
    if ndims(b) == 1
        (axes(A,2),)
    else
        (axes(A,2), axes(b, 2))
    end
end

function propagate_shape(::typeof(inv), A)
    @oops shp = shape(A)
    @assert ndims(A) == 2 && reverse(shp) == shp "Inv called on a non-square matrix"
    shp
end

function scalarize(arr::ArrayOp, idx)
    @assert length(arr.output_idx) == length(idx)

    axs = ranges(arr)

    iidx = collect(keys(axs))
    contracted = setdiff(iidx, arr.output_idx)

    axes = [get_extents(axs[c]) for c in contracted]
    summed = if isempty(contracted)
        arr.expr
    else
        mapreduce(arr.reduce, Iterators.product(axes...)) do idx
            replace_by_scalarizing(arr.expr, Dict(contracted .=> idx))
        end
    end

    dict = Dict(oi => (unwrap(i) isa Symbolic ? unwrap(i) : get_extents(axs[oi])[i])
                for (oi, i) in zip(arr.output_idx, idx) if unwrap(oi) isa Symbolic)

    replace_by_scalarizing(summed, dict)
end

scalarize(arr::Arr, idx) = wrap(scalarize(unwrap(arr),
                                          unwrap.(idx)))


eval_array_term(op, arr) = arr
eval_array_term(op::typeof(inv), arr) = inv(scalarize(wrap(arguments(arr)[1]))) #issue 653
eval_array_term(op::Arr) = wrap(eval_array_term(unwrap(op)))
eval_array_term(op) = eval_array_term(operation(op), op)

function scalarize(arr)
    if arr isa Arr || arr isa Symbolic{<:AbstractArray}
        if istree(arr)
            arr = eval_array_term(arr)
        end
        map(Iterators.product(axes(arr)...)) do i
            scalarize(arr[i...]) # Use arr[i...] here to trigger any getindex hooks
        end
    elseif istree(arr) && operation(arr) == getindex
        args = arguments(arr)
        scalarize(args[1], (args[2:end]...,))
    elseif arr isa Num
        wrap(scalarize(unwrap(arr)))
    elseif istree(arr) && symtype(arr) <: Number
        t = similarterm(arr, operation(arr), map(scalarize, arguments(arr)), symtype(arr), metadata=metadata(arr))
        istree(t) ? scalarize_op(operation(t), t) : t
    else
        arr
    end
end

@wrapped Base.isempty(x::AbstractArray) = shape(unwrap(x)) !== Unknown() && _iszero(length(x))
Base.collect(x::Arr) = scalarize(x)
Base.collect(x::SymArray) = scalarize(x)
isarraysymbolic(x) = unwrap(x) isa Symbolic && SymbolicUtils.symtype(unwrap(x)) <: AbstractArray

Base.convert(::Type{<:Array{<:Any, N}}, arr::Arr{<:Any, N}) where {N} = scalarize(arr)


### Stencils

struct ArrayMaker{T, AT<:AbstractArray} <: Symbolic{AT}
    shape
    sequence
    metadata
end

function ArrayMaker(a::ArrayLike; eltype=eltype(a))
    ArrayMaker{eltype}(size(a), Any[axes(a) => a])
end

function arraymaker(T, shape, views, seq...)
    ArrayMaker{T}(shape, [(views .=> seq)...], nothing)
end

istree(x::ArrayMaker) = true
operation(x::ArrayMaker) = arraymaker
arguments(x::ArrayMaker) = [eltype(x), shape(x), map(first, x.sequence), map(last, x.sequence)...]

shape(am::ArrayMaker) = am.shape

function ArrayMaker{T}(sz::NTuple{N, Integer}, seq::Array=[]; atype=Array, metadata=nothing) where {N,T}
    ArrayMaker{T, atype{T, N}}(map(x->1:x, sz), seq, metadata)
end

(::Type{ArrayMaker{T}})(i::Int...; atype=Array) where {T} = ArrayMaker{T}(i, atype=atype)

function Base.show(io::IO, ac::ArrayMaker)
    print(io, Expr(:call, :ArrayMaker, ac.shape,
                   Expr(:block, ac.sequence...)))
end

function get_indexers(ex)
    @assert ex.head == :ref
    arr = ex.args[1]
    args = map(((i,x),)->x == Symbol(":") ? :(1:lastindex($arr, $i)) : x, enumerate(ex.args[2:end]))
    replace_ends(arr, args)
end

function search_and_replace(expr, key, val)
    isequal(expr, key) && return val

    expr isa Expr ?
        Expr(expr.head, map(x->search_and_replace(x, key,val), expr.args)...) :
        expr
end

function replace_ends(arr, idx)
    [search_and_replace(ix, :end, :(lastindex($arr, $i)))
     for (i, ix) in enumerate(idx)]
end

macro setview!(definition, arrayop)
    setview(definition, arrayop, true)
end

macro setview(definition, arrayop)
    setview(definition, arrayop, false)
end

output_index_ranges(c::CartesianIndices) = c.indices
output_index_ranges(ix...) = ix

function setview(definition, arrayop, inplace)
    output_view = get_indexers(definition)
    output_ref = definition.args[1]

    function check_assignment(vw, op)
        try Base.Broadcast.broadcast_shape(map(length, vw), size(op))
        catch err
            if err isa DimensionMismatch
                throw(DimensionMismatch("setview did not work while assigning indices " *
                                        "$vw to $op. LHS has size $(map(length, vw)) "*
                                        "and RHS has size $(size(op)) " *
                                        "-- they need to be broadcastable."))
            else
                rethrow(err)
            end
        end
    end

    function push(inplace)
        if inplace
            function (am, vw, op)
                check_assignment(vw, op)
                # assert proper size match
                push!(am.sequence, vw => op)
                am
            end
        else
            function (am, vw, op)
                check_assignment(vw, op)
                if am isa ArrayMaker
                    typeof(am)(am.shape, vcat(am.sequence, vw => op))
                else
                    am = ArrayMaker(am)
                    push!(am.sequence, vw => op)
                    am
                end
            end
        end
    end
    quote
        $(push(inplace))($output_ref,
                         $output_index_ranges($(output_view...)), $unwrap($arrayop))
    end |> esc
end

macro makearray(definition, sequence)
    output_shape = get_indexers(definition)
    output_name = definition.args[1]

    seq = map(filter(x->!(x isa LineNumberNode), sequence.args)) do pair
        @assert pair.head == :call && pair.args[1] == :(=>)
        # TODO: make sure the same symbol is used for the lhs array
        :(@setview! $(pair.args[2]) $(pair.args[3]))
    end

    quote
        $output_name = $ArrayMaker{Real}(map(length, ($(output_shape...),)))
        $(seq...)
        $output_name = $wrap($output_name)
    end |> esc
end

function best_order(output_idx, ks, rs)
    unique!(filter(issym, vcat(reverse(output_idx)..., collect(ks))))
end

function _cat(x, xs...; dims)
    arrays = (x, xs...)
    if dims isa Integer
        sz = Base.cat_size_shape(Base.dims2cat(dims), arrays...)
        T = reduce(promote_type, eltype.(xs), init=eltype(x))
        newdim = cumsum(map(a->size(a, dims), arrays))
        start = 1
        A = ArrayMaker{T}(sz...)
        for (dim, array) in zip(newdim, arrays)
            idx = CartesianIndices(ntuple(n -> n==dims ?
                                          (start:dim) : (1:sz[n]), length(sz)))
            start = dim + 1

            @setview! A[idx] array
        end
        return A
    else
        error("Block diagonal concatenation not supported")
    end
end

# Base.cat(x::Arr, xs...; dims) = _cat(x, xs...; dims)
# Base.cat(x::AbstractArray, y::Arr, xs...; dims) = _cat(x, y, xs...; dims)

# vv uncomment these for a major release
# Base.vcat(x::Arr, xs::AbstractVecOrMat...) = cat(x, xs..., dims=1)
# Base.hcat(x::Arr, xs::AbstractVecOrMat...) = cat(x, xs..., dims=2)
# Base.vcat(x::AbstractVecOrMat, y::Arr, xs::AbstractVecOrMat...) = _cat(x, y, xs..., dims=1)
# Base.hcat(x::AbstractVecOrMat, y::Arr, xs::AbstractVecOrMat...) = _cat(x, y, xs..., dims=2)
# Base.vcat(x::Arr, y::Arr) = _cat(x, y, dims=1) # plug ambiguity
# Base.hcat(x::Arr, y::Arr) = _cat(x, y, dims=2)

function scalarize(x::ArrayMaker)
    T = eltype(x)
    A = Array{wrapper_type(T)}(undef, size(x))
    for (vw, arr) in x.sequence
        if any(x->x isa AbstractArray, vw)
            A[vw...] .= scalarize(arr)
        else
            A[vw...] = scalarize(arr)
        end
    end
    A
end

function scalarize(x::ArrayMaker, idx)
    for (vw, arr) in reverse(x.sequence) # last one wins
        if any(x->issym(x) || istree(x), idx)
            return term(getindex, x, idx...)
        end
        if all(in.(idx, vw))
            if symtype(arr) <: AbstractArray
                # Filter out non-array indices because the RHS will be one dim less
                el = [searchsortedfirst(v, i)
                      for (v, i) in zip(vw, idx) if v isa AbstractArray]
                return scalarize(arr[el...])
            else
                return arr
            end
        end
    end
    if !any(x->issym(x) || istree(x), idx) && all(in.(idx, axes(x)))
        throw(UndefRefError())
    end

    throw(BoundsError(x, idx))
end


### Codegen

function SymbolicUtils.Code.toexpr(x::ArrayOp, st)
    haskey(st.symbolify, x) && return st.symbolify[x]

    if istree(x.term)
        toexpr(x.term, st)
    else
        _array_toexpr(x, st)
    end
end

function SymbolicUtils.Code.toexpr(x::Arr, st)
    toexpr(unwrap(x), st)
end

function SymbolicUtils.Code.toexpr(x::ArrayMaker, st)
    _array_toexpr(x, st)
end

function _array_toexpr(x, st)
    outsym = Symbol("_out")
    N = length(shape(x))
    ex = :(let $outsym = zeros(Float64, map(length, ($(shape(x)...),)))
          $(inplace_expr(x, outsym))
          $outsym
      end) |> LiteralExpr
    toexpr(ex, st)
end

function inplace_expr(x, out_array, dict=nothing)
    x = unwrap(x)
    if symtype(x) <: Number
        :($out_array .= $x)
    else
        :($copy!($out_array, $x))
    end
end

function inplace_expr(x::ArrayMaker, out, dict=Dict())
    ex = []

    intermediates = Dict()
    for (i, (vw, op)) in enumerate(x.sequence)
        out′ = Symbol(out, "_", i)
        push!(ex, :($out′ = $view($out, $(vw...))))
        push!(ex, inplace_expr(unwrap(op), out′, intermediates))
    end

    Expr(:block, (:($sym = $ex) for (ex, sym) in  intermediates)..., ex...)
end

function inplace_expr(x::AbstractArray, out, intermediates=Dict())
    # TODO: extract more intermediates
    :(begin
          $([:($out[$(Tuple(idx)...)] = $(substitute(x, intermediates)[Tuple(idx)...])) for idx in eachindex(x)]...)
      end)
end

function inplace_builtin(term, outsym)
    isarr(n) = x->symtype(x) <: AbstractArray{<:Any, n}
    if istree(term) && operation(term) == (*) && length(arguments(term)) == 2
        A, B = arguments(term)
        isarr(2)(A) && (isarr(1)(B) || isarr(2)(B)) && return :($mul!($outsym, $A, $B))
    end
    return nothing
end

function find_inter(acc, expr)
    if !issym(expr) && symtype(expr) <: AbstractArray
        push!(acc, expr)
    elseif istree(expr)
        foreach(x -> find_inter(acc, x), arguments(expr))
    end
    acc
end

function get_inputs(x::ArrayOp)
    unique(find_inter([], x.expr))
end

function similar_arrayvar(ex, name)
    Sym{symtype(ex)}(name) #TODO: shape?
end

function reset_to_one(range)
    @assert step(range) == 1
    Base.OneTo(length(range))
end

function reset_sym(i)
    Sym{Int}(Symbol(nameof(i), "′"))
end

function inplace_expr(x::ArrayOp, outsym = :_out, intermediates = nothing)
    if x.term !== nothing
        ex = inplace_builtin(x.term, outsym)
        if ex !== nothing
            return ex
        end
    end

    rs = copy(ranges(x))

    inters = filter(!issym, get_inputs(x))
    intermediate_exprs = map(enumerate(inters)) do (i, ex)
        if !isnothing(intermediates)
            if haskey(intermediates, ex)
                return ex => intermediates[ex]
            else
                sym = similar_arrayvar(ex, Symbol(outsym, :_input_, i))
                intermediates[ex] = sym
                return ex => sym
            end
        else
            return ex => similar_arrayvar(ex, Symbol(outsym, :_input_, i))
        end
    end

    loops = best_order(x.output_idx, keys(rs), rs)

    expr = substitute(unwrap(x.expr), Dict(intermediate_exprs))

    out_idxs = map(reset_sym, x.output_idx)
    inner_expr = :($outsym[$(out_idxs...)] = $(x.reduce)($outsym[$(out_idxs...)], $(expr)))


    loops = foldl(reverse(loops), init=inner_expr) do acc, k
        if any(isequal(k), x.output_idx)
            :(for ($k, $(reset_sym(k))) in zip($(get_extents(rs[k])),
                                               reset_to_one($(get_extents(rs[k]))))
                  $acc
              end)
        else
            :(for $k in $(get_extents(rs[k]))
                  $acc
              end)
        end
    end

    if intermediates === nothing
        # output the intermediate generation
        :($(map(x->:($(x[2]) = $(x[1])), intermediate_exprs)...);
          $loops) |> SymbolicUtils.Code.LiteralExpr
    else
        SymbolicUtils.Code.LiteralExpr(loops)
    end
end


#=
"""
Find any inputs to ArrayOp that are ArrayMaker, and return
how to split all the inputs simultaneously so that the blocks
can now interact.
"""
function get_simultaneous_ranges(ex::ArrayOp)
    rs = ranges(ex)
    combine_together = []
    for (i, arrs) in rs
        together = unique(map(a->(a.A, a.dim), arrs))
        if length(together) > 1
            push!(combine_together, together)
        end
    end

    splits = map(combine_together) do group
        map(group) do a
            (A, dim) = a
            if A isa ArrayMaker
                sort(map(x->x[dim], map(first, A.sequence)), by=first)
            else
                [axes(A, dim)]
            end
        end
    end

    combined_splits = map(splits) do rs
        new_starts = sort!(unique!(reduce(vcat, map(x->first.(x), rs))))
        lst = maximum(map(maximum, map(x->last.(x), rs)))
        UnitRange.(new_starts, vcat((new_starts .- 1)[2:end], lst))
    end

    collected = Dict(A => Any[[1:size(A, dim)] for dim in 1:ndims(A)]
                     for A in unique(reduce(vcat, map(x->map(a->a.A, x),
                                                      collect(values(rs))))))

    for (dims, rs) in zip(combine_together, combined_splits)
        for d in dims
            collected[d[1]][d[2]] = rs
        end
    end
    collected
end
=#
