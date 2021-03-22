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

        ArrayShape(map(get, shp))
    end
end

function maybefoldl(f, g, xs, acc)
    for x in xs
        y = f(x)
        y === Unknown() && return Unknown()
        acc = g(acc, y)
    end
    return acc
end

function Base.adjoint(A::SymArray)
    N = nd(A)
    if N !== nothing && !(N in (1, 2))
        error("Can adjoint only a vector or a matrix")
    end

    arrterm(adjoint, A)

    shp = shape(A)
    if shp !== nothing
        shp = ArrayShape(ndims(shp) == 2 ? reverse(shp.axes) : (Base.OneTo(1), shp.axes...))
        A_t = similararraytype(symtype(A), N=2)
        setmetadata(Term{A_t}(adjoint, [A]),
                    ArrayShapeCtx,
                    shp
                    )
    else
        Term{symtype(A)}(adjoint, [A])
    end
end

import Base: *, \
function (*)(A::Symbolic{<:AbstractMatrix},
             b::Symbolic{<:AbstractVector})
    @syms i::Int k::Int
    if istree(A) &&
        operation(A) === adjoint &&
        nd(arguments(A)[1]) === 1
        # do this to fail if dim mismatch
        shape_propagate(TensorOp((i,), A[i,k] * b[k]))
        return Term{Real}(*, [A, b])
    end

    shp = shape_propagate(TensorOp((i,), A[i,k] * b[k]))
    setmetadata(Term{Vector}(*, [A, b]), ArrayShapeCtx, ArrayShape(map(get, shp)))
end

#=
# a .+ 1
#
# a --> Symarray with known size
#     a --> Symarray with size (m, n) (n, m)
# a --> Symarray with known dimension but no size
# a --> Sym{AbstractArray} without any shape info
=#
