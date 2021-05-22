using Symbolics, Test
using DomainSets
@variables x y

I1 = Integral( x , ClosedInterval(1, 5))
I2 = Integral( x , ClosedInterval(1, 5))
@test I1 == I2
