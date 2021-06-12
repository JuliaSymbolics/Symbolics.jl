using Symbolics, Test
using DomainSets
@variables x y

I1 = Integral( x , ClosedInterval(1, 5))
I2 = Integral( x , ClosedInterval(1, 5))
@test I1 == I2

@variables  v(..) , u(..) , x , y
D = Differential(x)
Dxx = Differential(x)^2
I = Integral(x, ClosedInterval(1 , v(x)))
eq = D((I(u(x,t)) )) ~ 0
O = expand_derivatives(eq.lhs)
println(O)
eq = Dxx((I(u(x,t)) )) ~ 0
O = expand_derivatives(eq.lhs)
println(O)
