using Symbolics
using Test

@variables t
D1 = Difference(t; dt=0.01)
D2 = Difference(t; dt=0.01)

@test D1 == D2
@test Base.isequal(D1, D2)
@test Base.hash(D1) == Base.hash(D2)
