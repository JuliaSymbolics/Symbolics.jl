module SymbolicsLuxCoreExt

using LuxCore, Symbolics

@register_array_symbolic LuxCore.stateless_apply(
    model::LuxCore.AbstractExplicitLayer, x::AbstractArray, ps::Union{NamedTuple, <:AbstractVector}) begin
    size = LuxCore.outputsize(model)
    eltype = Real
end

end
