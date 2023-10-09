using Symbolics
using LogExpFunctions

N = 10

@variables a, b, c, x[1:N]

_a = -0.2
_b = -1.0
_c = 2.0
_x = rand(N)

vals = Dict(a => _a, b => _b, c => _c, x => _x)

@test substitute(logaddexp(a, b), vals) ≈ logaddexp(_a, _b)
@test substitute(logaddexp(a, _b), vals) ≈ logaddexp(_a, _b)
@test substitute(logaddexp(_a, b), vals) ≈ logaddexp(_a, _b)
@test substitute(logsubexp(a, b), vals) ≈ logsubexp(_a, _b) 
@test substitute(logsubexp(a, _b), vals) ≈ logsubexp(_a, _b) 
@test substitute(logsubexp(_a, b), vals) ≈ logsubexp(_a, _b) 
@test substitute(log1mexp(a), vals) ≈ log1mexp(_a)
@test substitute(log1pexp(a), vals) ≈ log1pexp(_a)
@test substitute(logexpm1(c), vals) ≈ logexpm1(_c)
@test substitute(logmxp1(c), vals) ≈ logmxp1(_c)
@test substitute(logsumexp(x), vals) ≈ logsumexp(_x)