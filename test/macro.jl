import Symbolics
import Symbolics: getdefaultval
import SymbolicUtils: Term, symtype, FnType
using Test

@test !(@isdefined Num)
@variables t
Symbolics.@register fff(t)
@test isequal(fff(t), Symbolics.Num(Symbolics.Term{Real}(fff, [Symbolics.value(t)])))

@variables t=0 a=1 x[1:4]=2 y[1:4](t)=3 w[1:4] = 1:4 z[1:4](t) = 2:5 p[1:4](..)

@test getdefaultval(t) == 0
@test getdefaultval(a) == 1
@test_throws ErrorException getdefaultval(x)
@test getdefaultval(x[1]) == 2
@test getdefaultval(y[2]) == 3
@test getdefaultval(w[2]) == 2
@test getdefaultval(w[4]) == 4
@test getdefaultval(z[3]) == 4

@test p[1] isa Term
@test symtype(p[1]) <: FnType{Tuple, Real}
@test p[1](t) isa Num
