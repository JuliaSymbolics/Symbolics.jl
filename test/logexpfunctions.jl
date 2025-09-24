using Symbolics
using LogExpFunctions
using Test

N = 10

@variables a, b, c, x[1:N]

_a = -0.2
_b = -1.0
_c = 2.0
_x = rand(N)

vals = Dict(a => _a, b => _b, c => _c, x => _x)

@test substitute(logaddexp(a, b), vals; fold = Val(true)) ≈ Num(logaddexp(_a, _b)) atol=1e-15
@test substitute(logaddexp(a, _b), vals; fold = Val(true)) ≈ Num(logaddexp(_a, _b)) atol=1e-15
@test substitute(logaddexp(_a, b), vals; fold = Val(true)) ≈ Num(logaddexp(_a, _b)) atol=1e-15
@test substitute(logsubexp(a, b), vals; fold = Val(true)) ≈ Num(logsubexp(_a, _b) ) atol=1e-15
@test substitute(logsubexp(a, _b), vals; fold = Val(true)) ≈ Num(logsubexp(_a, _b) ) atol=1e-15
@test substitute(logsubexp(_a, b), vals; fold = Val(true)) ≈ Num(logsubexp(_a, _b) ) atol=1e-15
@test substitute(log1mexp(a), vals; fold = Val(true)) ≈ Num(log1mexp(_a)) atol=1e-15
@test substitute(log1pexp(a), vals; fold = Val(true)) ≈ Num(log1pexp(_a)) atol=1e-15
@test substitute(logexpm1(c), vals; fold = Val(true)) ≈ Num(logexpm1(_c)) atol=1e-15
@test substitute(logmxp1(c), vals; fold = Val(true)) ≈ Num(logmxp1(_c)) atol=1e-15
@test substitute(logsumexp(x), vals; fold = Val(true)) ≈ Num(logsumexp(_x)) atol=1e-15
