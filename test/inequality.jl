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

@test Symbolics.substitute(a ≲ 2, a => 1) == (1 ≲ 2)

@test (a ≳ 2) == Inequality(a, 2, Symbolics.geq)
@test (a * x ≳ b / z) == Inequality(a * x, b / z, Symbolics.geq)
@test (a ≳ u) == Inequality(a, u, Symbolics.geq)
@test (a ≳ sin(u)) == Inequality(a, sin(u), Symbolics.geq)
@test Symbolics.scalarize(v .≳ u) == [Inequality(v[1], u, Symbolics.geq), Inequality(v[2], u, Symbolics.geq), Inequality(v[3], u, Symbolics.geq)]
@test Symbolics.scalarize(v .≳ w .+ 3) == [Inequality(v[1], w[1] + 3, Symbolics.geq), Inequality(v[2], w[2] + 3, Symbolics.geq), Inequality(v[3], w[3] + 3, Symbolics.geq)]

@test Symbolics.canonical_form(a + b *c ≳ x + 2 * x) == (3x - a - b*c ≲ 0)

@test Symbolics.substitute(a ≳ 2, a => 1) == (1 ≳ 2)

@testset "Array operands" begin
    @variables p[1:3] q[1:3] M[1:2, 1:3] s
    for (rel, op) in ((≲, Symbolics.leq), (≳, Symbolics.geq))
        @test rel(p, q) == Inequality(p, q, op)
        @test rel(p, s) == Inequality(p, s, op)
        @test rel(s, p) == Inequality(s, p, op)
        @test rel(p, 1) == Inequality(p, 1, op)
        @test rel(v, u) == Inequality(v, u, op)
        @test rel([a, b], 1) == Inequality([a, b], 1, op)
        @test rel(M, zeros(2, 3)) == Inequality(M, zeros(2, 3), op)
        @test_throws ArgumentError rel(p', q)
        @test_throws ArgumentError rel(M, p)
        @test_throws ArgumentError rel([a, b], p)

        @test Symbolics.scalarize(rel(p, q)) == [Inequality(p[i], q[i], op) for i in 1:3]
        @test Symbolics.scalarize(rel(p, q)) == Symbolics.scalarize(broadcast(rel, p, q))
        @test Symbolics.scalarize(rel(p, s)) == [Inequality(p[i], s, op) for i in 1:3]
        @test Symbolics.scalarize(rel(s, p)) == [Inequality(s, p[i], op) for i in 1:3]
        @test Symbolics.scalarize(rel(M, s)) == [Inequality(M[i, j], s, op) for i in 1:2, j in 1:3]
    end

    reduction = sum(abs2, M .* p'; dims = 1) ≲ s
    @test reduction isa Inequality
    @test size(reduction.lhs) == (1, 3)

    @test isequal(Symbolics.scalarize(Symbolics.canonical_form(p ≲ s).lhs), [p[i] - s for i in 1:3])
    @test isequal(Symbolics.scalarize(Symbolics.canonical_form(p ≳ s).lhs), [s - p[i] for i in 1:3])
    @test isequal(Symbolics.scalarize(Symbolics.canonical_form(p ≲ q).lhs), [p[i] - q[i] for i in 1:3])
    @test Symbolics.canonical_form(p ≳ q).relational_op == Symbolics.leq

    @test Symbolics.substitute(p ≲ s, Dict(s => 2)) == (p ≲ 2)

    @test Symbolics.evaluate(p ≲ q, Dict(p => [1, 2, 3], q => [2, 3, 4]))
    @test !Symbolics.evaluate(p ≲ q, Dict(p => [1, 2, 3], q => [2, 0, 5]))
    @test Symbolics.evaluate(p ≳ q, Dict(p => [2, 3, 4], q => [1, 2, 3]))
    @test !Symbolics.evaluate(p ≳ q, Dict(p => [2, 3, 4], q => [1, 5, 3]))
    @test Symbolics.evaluate(p ≲ s, Dict(p => [1, 2, 3], s => 4))
    @test !Symbolics.evaluate(p ≲ s, Dict(p => [1, 5, 3], s => 4))
end
