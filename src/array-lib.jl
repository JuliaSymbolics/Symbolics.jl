#################### BROADCAST ################
#
using Base.Broadcast

Base.broadcastable(s::Symbolic{<:AbstractArray}) = s
struct SymBroadcast <: Broadcast.BroadcastStyle end
Broadcast.BroadcastStyle(::Type{<:Symbolic{<:AbstractArray}}) = SymBroadcast()
Broadcast.result_style(::SymBroadcast) = SymBroadcast()
Broadcast.BroadcastStyle(::SymBroadcast, ::Broadcast.BroadcastStyle) = SymBroadcast()

function Broadcast.materialize(bc::Broadcast.Broadcasted{SymBroadcast})
    # Do the thing here
    arrterm(broadcast, bc.f, bc.args...)
end

function propagate_ndims(::typeof(broadcast), f, args...)
    maybefoldl(getndims, max, args, 0)
end

function propagate_shape(::typeof(broadcast), f, args...)
    ndim = propagate_ndims(broadcast, f, args...)
    maybe(ndim) do ndim
        subscripts = [Sym{Int}(Symbol("i_$i")) for i in 1:ndim]
        args′ = map(args) do x
            if ndims(x) != 0
                subs = map(i-> isone(size(x, i)) ? 1 : subscripts[i], 1:ndims(x))
                term(getindex, x, subs...)
            else
                x
            end
        end

        shp = shape_propagate(TensorOp((subscripts...,),
                                       term(+, args′...)))

        map(get, shp)
    end
end
# propagate_atype, propagate_eltype
#################### TRANSPOSE ################
#
function Base.adjoint(A::SymArray)
    N = getndims(A)
    if N !== Unknown() && !(N in (1, 2))
        error("Can adjoint only a vector or a matrix")
    end

    arrterm(adjoint, A)
end

propagate_ndims(::typeof(adjoint), A) = 2
function propagate_shape(::typeof(adjoint), A)
    @oops shp = shape(A)
    length(shp) == 2 ? reverse(shp) : (Base.OneTo(1), shp...)
end

#################### MATMUL ################
#
import Base: *, \
function (*)(A::Symbolic{<:AbstractMatrix},
             b::Union{Symbolic{<:AbstractMatrix},
                      Symbolic{<:AbstractVector}})
    if istree(A) &&
        operation(A) === adjoint &&
        getndims(arguments(A)[1]) === 1
        # do this to fail if dim mismatch
        T = Base.promote_op(*, eltype(A), eltype(b))
        T = T == Any ? Number : T
        Term{T}(*, [A, b])
    else
        arrterm(*, A, b)
    end
end

propagate_ndims(::typeof(*), A, B) = getndims(B)
function propagate_shape(::typeof(*), A, b::Symbolic{<:AbstractVector})
    @syms i::Int k::Int
    shp = shape_propagate(TensorOp((i,), A[i,k] * b[k]))
    map(get, shp)
end

function propagate_shape(::typeof(*), A, B::Symbolic{<:AbstractMatrix})
    @syms i::Int j::Int k::Int
    shp = shape_propagate(TensorOp((i,j), A[i,k] * B[k,j]))
    map(get, shp)
end

#################### MAP-REDUCE ################
#

Base.map(f, x::SymArray) = arrterm(map, f, x)
Base.map(f, x::SymArray, xs...) = arrterm(map, f, x, xs...)
Base.map(f, x, y::SymArray, z...) = arrterm(map, f, x, y, z...)
Base.map(f, x, y, z::SymArray, w...) = arrterm(map, f, x, y, z, w...)
function propagate_eltype(::typeof(map), f, xs...)
    if any(x->x isa Unknown, geteltype.(xs))
        return Unknown()
    end
    Base.return_types(f, eltype.(xs))[1]
end

# TODO: find one that has a shape, not just the first one.
propagate_shape(::typeof(map), f, x, xs...) = shape(x)


@inline _mapreduce(f, x, dims, kw) = mapreduce(f, x; dims=dims, kw...)

function Base.mapreduce(f, g, x::SymArray; dims=:, kw...)
    if dims === (:)
        return Term{Number}(_mapreduce, [f, g, x, dims, (kw...,)])
    end

    arrterm(_mapreduce, f, g, x, dims, kw)
end

function propagate_shape(::typeof(_mapreduce), f, g, x, dims, kw)
    @oops shp = shape(x)

    map(enumerate(shp)) do (i, dim)
        i in dims ? Base.OneTo(1) : dim
    end
end

for f in [sum, prod, any, all]
    _f = Symbol("_", nameof(f))
    T = f in (any, all) ? Bool : Number

    @eval function $_f(x, dims, kw)
        $f(x; dims=dims, kw...)
    end

    @eval function (::$(typeof(f)))(x::SymArray; dims=:, kw...)
        if dims === (:)
            return Term{$T}($_f, [x, dims, (kw...,)])
        end

        arrterm($_f, x, dims, (kw...,))
    end

    @eval function propagate_shape(::typeof($_f), x, dims, kw)
        @oops shp = shape(x)

        (map(enumerate(shp)) do (i, dim)
            i in dims ? Base.OneTo(1) : dim
        end...,)
    end
end


#=
#
#Bounds check
#
julia> @syms x[1:10]
(x,)

julia> x[1:11]
getindex(x, 1:11)
=#
