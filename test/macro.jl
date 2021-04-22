import ModelingToolkit
using Test

@test !(@isdefined Num)
ModelingToolkit.@variables t
ModelingToolkit.@register fff(t)
@test isequal(fff(t), ModelingToolkit.Num(ModelingToolkit.Term{Real}(fff, [ModelingToolkit.value(t)])))
