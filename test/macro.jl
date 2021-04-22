import Symbolics
using Test

@test !(@isdefined Num)
Symbolics.@variables t
Symbolics.@register fff(t)
@test isequal(fff(t), Symbolics.Num(Symbolics.Term{Real}(fff, [Symbolics.value(t)])))
