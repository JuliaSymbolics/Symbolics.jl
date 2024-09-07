module SymbolicsLuxExt

using Lux
using Symbolics
using Lux.LuxCore
using Symbolics.SymbolicUtils

function Lux.NilSizePropagation.recursively_nillify(x::SymbolicUtils.BasicSymbolic{<:Vector{<:Real}})
    Lux.NilSizePropagation.recursively_nillify(Symbolics.wrap(x))
end

@register_array_symbolic LuxCore.stateless_apply(
    model::LuxCore.AbstractLuxLayer, x::AbstractArray, ps::Union{NamedTuple, <:AbstractVector}) begin
    size = LuxCore.outputsize(model, x, LuxCore.Random.default_rng())
    eltype = Real
end

end
