using Symbolics
using Symbolics: hasdiff
using Test

@variables t x
D1 = Difference(t; dt = 0.01)
D2 = Difference(t; dt = 0.01)

@test D1 == D2
@test Base.isequal(D1, D2)
@test Base.hash(D1) == Base.hash(D2)

@test D1(x) isa Num

@test isequal((D1^2)(x), D1(D1(x)))

# hasdiff
@test hasdiff(D1)
@test hasdiff(D1(x) ~ x)
@test hasdiff(x ~ D1(x))
@test !hasdiff(t ~ x)
@test hasdiff((D1^2)(x) ~ x)
@test hasdiff(D1(x) ~ D2(x))
