using Symbolics, Test
import Symbolics: VarDomainPairing
import DomainSets: Interval, Ball, infimum, supremum, radius, center

domain = Interval(0, 1)
@test infimum(domain) == 0
@test supremum(domain) == 1

ball = Ball()
@test radius(ball) == 1.0
@test center(ball) == [0.0, 0.0, 0.0]

ball = Ball(2.5, [1, 2])
@test radius(ball) == 2.5
@test center(ball) == [1, 2]

@syms t x
t in domain
domains = [t ∈ Interval(0.0, 1.0),
    x ∈ Interval(0.0, 1.0)]
@test domains[1] isa VarDomainPairing
@test domains[2] isa VarDomainPairing

@syms y z
@test ((x, y) ∈ Ball(2.0, [0, 0])) isa VarDomainPairing
@test ((x, y, z) ∈ Ball(1.5, [1, 2, 3])) isa VarDomainPairing

# Test tuple gets converted to Interval
var_domain_pair = t ∈ (0, 1)
@test var_domain_pair isa VarDomainPairing
@test var_domain_pair.domain isa Interval

# Other types
t = Symbolics.Num(:t)
@assert (t ∈ domain) isa VarDomainPairing
t = Symbolics.variable(:t)
@assert (t ∈ domain) isa VarDomainPairing
