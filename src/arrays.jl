using SymbolicUtils
using StaticArraysCore
import Base: eltype, length, ndims, size, axes, eachindex

### Wrapper type for dispatch

@symbolic_wrap struct Arr{T,N} <: AbstractArray{T, N}
    value::BasicSymbolic{VartypeT}

    function Arr{T, N}(ex) where {T, N}
        if is_wrapper_type(T)
            @assert symtype(ex) <: AbstractArray{<:wraps_type(T), N}
        else
            @assert symtype(ex) <: AbstractArray{T, N}
        end
        new{T, N}(Const{VartypeT}(ex))
    end
end

Base.hash(x::Arr, u::UInt) = hash(unwrap(x), u)
Base.isequal(a::Arr, b::Arr) = isequal(unwrap(a), unwrap(b))
Base.isequal(a::Arr, b::BasicSymbolic) = isequal(unwrap(a), b)
Base.isequal(a::BasicSymbolic, b::Arr) = isequal(b, a)

function Arr(x)
    A = symtype(x)
    @assert A <: AbstractArray
    Arr{maybewrap(eltype(A)), ndims(A)}(x)
end

SymbolicUtils.unwrap(x::Arr) = x.value
SymbolicUtils.symtype(x::Arr) = symtype(unwrap(x))

maybewrap(T) = has_symwrapper(T) ? wrapper_type(T) : T
# These methods allow @wrapped methods to be more specific and not overwrite
# each other when defined both for matrix and vector
wrapper_type(::Type{<:AbstractMatrix}) = Arr{<:Any, 2}
wrapper_type(::Type{<:AbstractMatrix{T}}) where {T} = Arr{maybewrap(T), 2}

wrapper_type(::Type{<:AbstractVector}) = Arr{<:Any, 1}
wrapper_type(::Type{<:AbstractVector{T}}) where {T} = Arr{maybewrap(T), 1}

function Base.show(io::IO, arr::Arr)
    x = unwrap(arr)
    iscall(x) && print(io, "(")
    print(io, unwrap(arr))
    iscall(x) && print(io, ")")
    print(io, "[")
    if shape(x) isa SymbolicUtils.Unknown
        print(io, shape(x))
    else
        print(io, join(string.(axes(arr)), ","))
    end
    print(io, "]")
end
Base.show(io::IO, ::MIME"text/plain", arr::Arr) = show(io, arr)

################# Base array functions
Base.length(A::Arr) = length(unwrap(A))
Base.size(A::Arr) = size(unwrap(A))
Base.axes(A::Arr) = axes(unwrap(A))
Base.eachindex(A::Arr) = eachindex(unwrap(A))

function SymbolicUtils.search_variables!(buffer, expr::Arr; kw...)
    SymbolicUtils.search_variables!(buffer, unwrap(expr); kw...)
end

# cannot use `@wrapped` since it will define `\(::BasicSymbolic, ::BasicSymbolic)`
# and because `\(::Arr, ::BasicSymbolic)` will be ambiguous.
for (T1, T2) in [
    (Arr{<:Any, 2}, Arr{<:Any, 1}),
    (Arr{<:Any, 2}, Arr{<:Any, 2}),
    (Arr{<:Any, 2}, BasicSymbolic{SymReal}),
    (Arr{<:Any, 2}, BasicSymbolic{SafeReal}),
    (Arr{<:Any, 2}, BasicSymbolic{TreeReal}),
    (BasicSymbolic{SymReal}, Arr{<:Any, 1}),
    (BasicSymbolic{SafeReal}, Arr{<:Any, 1}),
    (BasicSymbolic{TreeReal}, Arr{<:Any, 1}),
    (BasicSymbolic{SymReal}, Arr{<:Any, 2}),
    (BasicSymbolic{SafeReal}, Arr{<:Any, 2}),
    (BasicSymbolic{TreeReal}, Arr{<:Any, 2}),
]
    @eval function Base.:(\)(A::$T1, b::$T2)
        unwrap(A) \ unwrap(b)
    end
end

Base.inv(A::Arr{<:Any, 2}) = wrap(inv(unwrap(A)))
LinearAlgebra.det(A::Arr{<:Any, 2}) = wrap(det(unwrap(A)))
LinearAlgebra.adjoint(A::Arr{<:Any, 2}) = wrap(adjoint(unwrap(A)))
LinearAlgebra.adjoint(A::Arr{<:Any, 1}) = wrap(adjoint(unwrap(A)))

SymbolicUtils.scalarize(x::Arr) = SymbolicUtils.scalarize(unwrap(x))

Base.isempty(x::Arr) = isempty(unwrap(x))
Base.collect(x::Arr) = wrap.(collect(unwrap(x)))
isarraysymbolic(x) = false
# this should be validated in the constructor
isarraysymbolic(x::Arr) = true
isarraysymbolic(x::BasicSymbolic) = symtype(x) <: AbstractArray

Base.convert(::Type{<:Array{<:Any, N}}, arr::Arr{<:Any, N}) where {N} = scalarize(arr)

function SymbolicUtils.Code.toexpr(x::Arr, st)
    toexpr(unwrap(x), st)
end

"""
    $(TYPEDSIGNATURES)

Utility to represent `dst .= src` as a function call.
"""
function broadcast_assign!(dst, src)
    dst .= src
end

function inplace_expr(x, out_array, intermediates = nothing)
    x = unwrap(x)
    if SymbolicUtils.isarrayop(x)
        return x
    elseif symtype(x) <: Number
        return term(broadcast_assign!, out_array, x)
    else
        return term(copy!, out_array, x)
    end
end

function inplace_expr(x::AbstractArray, out, intermediates = nothing)
    expr = SetArray(false, out, x, true)
    if intermediates !== nothing
        expr = Let(map(k -> Assignment(intermediates[k], k), collect(keys(intermediates))), expr, false)
    end
    return expr
end

function inplace_builtin(term, outsym)
    isarr(n) = x->symtype(x) <: AbstractArray{<:Any, n}
    if iscall(term) && operation(term) == (*) && length(arguments(term)) == 2
        A, B = arguments(term)
        isarr(2)(A) && (isarr(1)(B) || isarr(2)(B)) && return term(mul!, outsym, A, B)
    end
    return nothing
end

hasnode(r::Function, y::Arr) = _hasnode(r, y)
hasnode(r::Union{Num, BasicSymbolic, Arr}, y::Arr) = occursin(unwrap(r), unwrap(y))

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
