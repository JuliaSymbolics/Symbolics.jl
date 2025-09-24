module SymbolicsLuxExt

using Lux
using Symbolics
using Lux.LuxCore
using Lux.Random: AbstractRNG, default_rng
using Symbolics.SymbolicUtils

@static if isdefined(Lux.NilSizePropagation, :recursively_nillify)
    function Lux.NilSizePropagation.recursively_nillify(x::SymbolicUtils.BasicSymbolic)
        @assert SymbolicUtils.symtype(x) <: Vector{<:Real}
        Lux.NilSizePropagation.recursively_nillify(Symbolics.wrap(x))
    end
end

function LuxCore.outputsize(model::SymbolicUtils.BasicSymbolic, x::Symbolics.Arr, rng::AbstractRNG)
    @assert SymbolicUtils.symtype(model) <: LuxCore.AbstractLuxLayer
    concrete_model = if SymbolicUtils.isconst(model)
        SymbolicUtils.unwrap_const(model)
    else
        Symbolics.getdefaultval(model)
    end
    LuxCore.outputsize(concrete_model, x, rng)
end

@register_array_symbolic LuxCore.stateless_apply(
    model::LuxCore.AbstractLuxLayer, x::AbstractArray, ps::Union{NamedTuple, <:AbstractVector}) begin
    size = LuxCore.outputsize(model, Symbolics.wrap(x), default_rng())
    eltype = Real
end

end
