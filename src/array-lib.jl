using Base.Broadcast

Base.broadcastable(s::Symbolic{<:AbstractArray}) = s
struct SymBroadcast <: Broadcast.BroadcastStyle end
Broadcast.BroadcastStyle(::Type{<:Symbolic{<:AbstractArray}}) = SymBroadcast()
Broadcast.result_style(::SymBroadcast) = SymBroadcast()
Broadcast.BroadcastStyle(::SymBroadcast, ::Broadcast.BroadcastStyle) = SymBroadcast()


function bestarraytype(A)
    @maybe T=elt(A) begin
        @maybe N=nd(A) slicetype(A){T,N}
        return slicetype(A){T}
    end

    @maybe N=nd(A) return slicetype(A){N}

    return slicetype(A)
end

function Broadcast.materialize(bc::Broadcast.Broadcasted{SymBroadcast})
    # Do the thing here
    ndim = maybefoldl(nd, max, bc.args, 0)
    if ndim === nothing
        return Term{AbstractArray}(broadcast, [bc.f, bc.args...])
    else
        subscripts = [Sym{Int}(Symbol("i_$i")) for i in 1:ndim]
        args = map(bc.args) do x
            if ndims(x) != 0
                subs = map(i-> isone(size(x, i)) ? 1 : subscripts[i], 1:ndims(x))
                term(getindex, x, subs...)
            else
                x
            end
        end

        shp = shape_propagate(TensorOp((subscripts...,),
                                       term(+, args...)))

        setmetadata(Term{AbstractArray}(broadcast,
                                        [bc.f, bc.args...]),
                    ArrayShapeCtx,
                    ArrayShape(map(get, shp)))
    end
end

function maybefoldl(f, g, xs, acc)
    for x in xs
        @maybe y=f(x) begin
            acc = g(acc, y)
            continue
        end
        return nothing
    end
    return acc
end


function Base.adjoint(A::SymArray)
    N = nd(A)
    if N !== nothing && !(N in (1, 2))
        error("Can adjoint only a vector or a matrix")
    end

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

#=
# a .+ 1
#
# a --> Symarray with known size
#     a --> Symarray with size (m, n) (n, m)
# a --> Symarray with known dimension but no size
# a --> Sym{AbstractArray} without any shape info
=#
