using Symbolics, Test
import DomainSets: Interval, Ball, infimum, supremum, radius, center

domain = Interval(0, 1)
@test infimum(domain) == 0
@test supremum(domain) == 1

ball = Ball()
@test radius(ball) == 1.0
@test center(ball) == [0.0,0.0,0.0]

ball = Ball(2.5, [1,2])
@test radius(ball) == 2.5
@test center(ball) == [1,2]

@syms t x
t in domain
domains = [t ∈ Interval(0.0,1.0),
           x ∈ Interval(0.0,1.0)]

# Test tuple gets converted to Interval
var_domain_pair = t ∈ (0,1)
@test var_domain_pair.domain isa Interval

@syms y z
(x,y) ∈ Ball(2.0, [0,0])
(x,y,z) ∈ Ball(1.5, [1,2,3])