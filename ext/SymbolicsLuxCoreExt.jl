module SymbolicsLuxCoreExt

using LuxCore, Symbolics

@register_array_symbolic LuxCore.partial_apply(
    model::LuxCore.AbstractExplicitContainerLayer, x::AbstractArray, ps::NamedTuple, st::NamedTuple) begin
    size = ((model[end].out_dims),)
    eltype = Real
end

@register_array_symbolic LuxCore.partial_apply(
    model::LuxCore.AbstractExplicitLayer, x::AbstractArray, ps::NamedTuple, st::NamedTuple) begin
    size = ((model.out_dims),)
    eltype = Real
end

end
