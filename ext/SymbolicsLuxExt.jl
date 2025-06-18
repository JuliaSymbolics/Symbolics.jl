module SymbolicsLuxExt

using Lux
using Symbolics
using Lux.LuxCore
using Lux.Random: AbstractRNG, default_rng
using Symbolics.SymbolicUtils

@static if isdefined(Lux.NilSizePropagation, :recursively_nillify)
    function Lux.NilSizePropagation.recursively_nillify(x::SymbolicUtils.BasicSymbolic{<:Vector{<:Real}})
        Lux.NilSizePropagation.recursively_nillify(Symbolics.wrap(x))
    end
end

function LuxCore.outputsize(model::SymbolicUtils.BasicSymbolic{<:LuxCore.AbstractLuxLayer}, x::Symbolics.Arr, rng::AbstractRNG)
    LuxCore.outputsize(Symbolics.getdefaultval(model), x, rng)
end

@register_array_symbolic LuxCore.stateless_apply(
    model::LuxCore.AbstractLuxLayer, x::AbstractArray, ps::Union{NamedTuple, <:AbstractVector}) begin
    size = LuxCore.outputsize(model, Symbolics.wrap(x), default_rng())
    eltype = Real
end

end
