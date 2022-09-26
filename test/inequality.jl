using Symbolics
using Test

@variables t a b c
@variables x y z
@variables u(t) v(t)[1:3] w(t)[1:3]

@test (a ≲ 2) == Inequality(a, 2, Symbolics.leq)
@test (a * x ≲ b / z) == Inequality(a * x, b / z, Symbolics.leq)
@test (a ≲ u) == Inequality(a, u, Symbolics.leq)
@test (a ≲ sin(u)) == Inequality(a, sin(u), Symbolics.leq)
@test Symbolics.scalarize(v .≲ u) == [Inequality(v[1], u, Symbolics.leq), Inequality(v[2], u, Symbolics.leq), Inequality(v[3], u, Symbolics.leq)]
@test Symbolics.scalarize(v .≲ w .+ 3) == [Inequality(v[1], w[1] + 3, Symbolics.leq), Inequality(v[2], w[2] + 3, Symbolics.leq), Inequality(v[3], w[3] + 3, Symbolics.leq)]

@test Symbolics.canonical_form(a + b *c ≲ x + 2 * x) == (a + b*c - 3x ≲ 0)

@test (a ≳ 2) == Inequality(a, 2, Symbolics.geq)
@test (a * x ≳ b / z) == Inequality(a * x, b / z, Symbolics.geq)
@test (a ≳ u) == Inequality(a, u, Symbolics.geq)
@test (a ≳ sin(u)) == Inequality(a, sin(u), Symbolics.geq)
@test Symbolics.scalarize(v .≳ u) == [Inequality(v[1], u, Symbolics.geq), Inequality(v[2], u, Symbolics.geq), Inequality(v[3], u, Symbolics.geq)]
@test Symbolics.scalarize(v .≳ w .+ 3) == [Inequality(v[1], w[1] + 3, Symbolics.geq), Inequality(v[2], w[2] + 3, Symbolics.geq), Inequality(v[3], w[3] + 3, Symbolics.geq)]

@test Symbolics.canonical_form(a + b *c ≳ x + 2 * x) == (3x - a - b*c ≲ 0)
