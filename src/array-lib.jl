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
    N = nd(A)
    if N !== nothing && !(N in (1, 2))
        error("Can adjoint only a vector or a matrix")
    end

    arrterm(adjoint, A)
end

propagate_ndims(::typeof(adjoint), A) = 2
function propagate_shape(::typeof(adjoint), A)
    shp = shape(A)
    shp === Unknown() && return Unknown()
    ndims(shp) == 2 ? reverse(shp) : (Base.OneTo(1), shp...)
end

#################### MATVEC ################
#
import Base: *, \
function (*)(A::Symbolic{<:AbstractMatrix},
             b::Symbolic{<:AbstractVector})
    arrterm(*, A, b)
end

propagate_ndims(::typeof(*), A, B) = getndims(B)
function propagate_shape(::typeof(*), A, B::Symbolic{<:AbstractVector})
    @syms i::Int k::Int
    if istree(A) &&
        operation(A) === adjoint &&
        nd(arguments(A)[1]) === 1
        # do this to fail if dim mismatch
        map(get, shape_propagate(TensorOp((i,), A[i,k] * b[k])))
    else
        shp = shape_propagate(TensorOp((i,), A[i,k] * b[k]))
        map(get, shp)
    end
end

function (*)(A::Symbolic{<:AbstractMatrix},
             B::Symbolic{<:AbstractMatrix})
    @syms i::Int j::Int k::Int
end

#=
# a .+ 1
#
# a --> Symarray with known size
#     a --> Symarray with size (m, n) (n, m)
# a --> Symarray with known dimension but no size
# a --> Sym{AbstractArray} without any shape info
=#
