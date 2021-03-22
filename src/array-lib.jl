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
                term(getindex, x, subscripts[1:ndims(x)]...)
            else
                x
            end
        end
        shp = shape_propagate(TensorOp((subscripts...,), term(+, args...)))
        setmetadata(Term{AbstractArray}(broadcast, [bc.f, bc.args...]), ArrayShapeCtx, ArrayShape(map(get, shp)))
    end
end

function promote_symtype(::typeof(broadcast),
                         f, args...)
    D = count(x->x <: Number, idx)
    @maybe T=elt(A) begin
        @maybe N=nd(A) return N-D == 0 ? T : slicetype(A){T,N-D}
        return slicetype(A){T}
    end

    @maybe N=nd(A) return N-D == 0 ? T : slicetype(A){T, N-D} where T

    return slicetype(A)
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

#=
# a .+ 1
#
# a --> Symarray with known size
#     a --> Symarray with size (m, n) (n, m)
# a --> Symarray with known dimension but no size
# a --> Sym{AbstractArray} without any shape info
=#
