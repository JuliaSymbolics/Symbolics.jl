using Symbolics
using Test

@variables t σ ρ β
@variables x y z
D = Difference(t; dt=0.01)

@test D == D
@test Base.isequal(D, D)
@test Base.hash(D) == Base.hash(D)
