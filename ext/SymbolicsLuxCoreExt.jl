module SymbolicsLuxCoreExt

using LuxCore, Symbolics

@register_array_symbolic LuxCore.stateless_apply(
    model::LuxCore.AbstractLuxLayer, x::AbstractArray, ps::Union{NamedTuple, <:AbstractVector}) begin
    size = LuxCore.outputsize(model, x, LuxCore.Lux.Random())
    eltype = Real
end

end
