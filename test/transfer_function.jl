using Symbolics
using Symbolics: value
using LinearAlgebra
using Random
using Test

iszero_expr(x) = isequal(value(simplify_fractions(x)), 0)
reference(A, B, C, D, ω) = C * ((im * ω * I - A) \ B) + D

@variables s R L C T m c k

# `(A, B, C, D)` in the form `ModelingToolkit.linearize_symbolic` returns them
@testset "first-order lag" begin
    G = transfer_function(fill(-1 / T, 1, 1), fill(1 / T, 1, 1), fill(Num(1), 1, 1), fill(Num(0), 1, 1), s)
    @test G isa Matrix{Num}
    @test size(G) == (1, 1)
    @test isequal(G[1], 1 / (1 + T * s))
end

@testset "series RLC" begin
    # states: capacitor voltage and inductor current
    A = [0 1 / C; -1 / L -R / L]
    B = reshape([0, 1 / L], 2, 1)
    G = transfer_function(A, B, Num[1 0], fill(Num(0), 1, 1), s)
    @test isequal(G[1], 1 / expand(L * C * s^2 + R * C * s + 1))
    # both states as outputs: the current is `C s` times the voltage
    G = transfer_function(A, B, Num[1 0; 0 1], zeros(Num, 2, 1), s)
    @test iszero_expr(G[2] - C * s * G[1])
end

@testset "mass-spring-damper" begin
    G = transfer_function([0 1; -k / m -c / m], [0, 1 / m], [1, 0], 0, s)
    @test isequal(G[1], 1 / expand(m * s^2 + c * s + k))
end

@testset "direct feedthrough" begin
    # RC high-pass
    τ = R * C
    G = transfer_function(fill(-1 / τ, 1, 1), [1 / τ], [-1], 1, s)
    @test iszero_expr(G[1] - τ * s / (τ * s + 1))
    # symbolic `D` with a denominator
    G = transfer_function(fill(-1 / T, 1, 1), [1 / T], [1], 1 / R, s)
    @test iszero_expr(G[1] - (R + T * s + 1) / (R * (T * s + 1)))
end

@testset "denominators are shared across a row" begin
    G = transfer_function([-1 / (R * C) 1 / C; 0 -1 / T], [0, 1 / T], [1, 0], 0, s)
    @test iszero_expr(G[1] - R / ((R * C * s + 1) * (T * s + 1)))
    num, den = Symbolics.arguments(Symbolics.unwrap(G[1]))
    @test Symbolics.degree(Num(den), C) == 1
end

@testset "sums of fractions in an entry" begin
    A = fill(-1 / (R * C) - 1 / (L * C), 1, 1)
    G = transfer_function(A, [1], [1], 0, s)
    @test iszero_expr(G[1] - 1 / (s + 1 / (R * C) + 1 / (L * C)))
end

@testset "cancellation" begin
    A = [-R 0; 0 -L]
    G = transfer_function(A, [1, 0], [1, 0], 0, s)
    @test iszero_expr(G[1] - 1 / (s + R))
    G = transfer_function(A, [1, 0], [1, 0], 0, s; simplify = true)
    @test isequal(G[1], simplify_fractions(1 / (s + R)))
    # decoupled input-output pairs are exactly zero
    G = transfer_function(A, [1 0; 0 1], [1 0; 0 1], zeros(2, 2), s)
    @test isequal(value(G[1, 2]), 0)
    @test isequal(value(G[2, 1]), 0)
end

@testset "symbolic array" begin
    @variables a[1:2, 1:2]
    G = transfer_function(a, [1, 0], [1, 0], 0, s)
    hand = (s - a[2, 2]) / ((s - a[1, 1]) * (s - a[2, 2]) - a[1, 2] * a[2, 1])
    @test iszero_expr(G[1] - hand)
    @test_throws ArgumentError frequency_response(G, s, [1.0])
end

@testset "6x6 sparse symbolic" begin
    @variables p[1:16]
    A = zeros(Num, 6, 6)
    for i in 1:6
        A[i, i] = -p[i]
        i < 6 || continue
        A[i, i + 1] = p[6 + i]
        A[i + 1, i] = -p[11 + i]
    end
    B = [1; zeros(Int, 5)]
    G = transfer_function(A, B, reverse(B), 0, s)
    vals = Dict(p[i] => 0.3 + i / 7 for i in 1:16)
    Anum = Float64.(value.(substitute.(A, Ref(vals); fold = Val(true))))
    ωs = [0.1, 1.0, 3.0]
    H = frequency_response(substitute.(G, Ref(vals)), s, ωs)
    for (i, ω) in enumerate(ωs)
        @test H[1, 1, i] ≈ reference(Anum, B, reverse(B)', 0, ω)
    end
end

@testset "random numeric systems" begin
    rng = MersenneTwister(1)
    ωs = range(0.0, 10.0, length = 7)
    for (nx, nu, ny) in ((1, 1, 1), (2, 3, 1), (3, 1, 2), (4, 2, 2), (6, 2, 3))
        A, B = randn(rng, nx, nx), randn(rng, nx, nu)
        C, D = randn(rng, ny, nx), randn(rng, ny, nu)
        H = frequency_response(transfer_function(A, B, C, D, s), s, ωs)
        @test H isa Array{ComplexF64, 3}
        @test size(H) == (ny, nu, length(ωs))
        for (i, ω) in enumerate(ωs)
            @test H[:, :, i] ≈ reference(A, B, C, D, ω) rtol = 1.0e-8
        end
    end
end

@testset "integer grid and complex-valued functions" begin
    G = transfer_function(fill(-1, 1, 1), [1], [1], 0, s)
    @test frequency_response(G, s, 0:2) ≈ reshape(1 ./ (im .* (0:2) .+ 1), 1, 1, 3)
    # a delay is not a state-space system, but the matrix can be evaluated
    H = frequency_response(reshape([exp(-s) * sqrt(s)], 1, 1), s, [2.0])
    @test H[1] ≈ exp(-2im) * sqrt(2im)
end

@testset "no states" begin
    D = Num[1 R; 0 L]
    G = transfer_function(zeros(0, 0), zeros(0, 2), zeros(2, 0), D, s)
    @test isequal(G, D)
    H = frequency_response(transfer_function(zeros(0, 0), zeros(0, 1), zeros(1, 0), 2, s), s, [0.0, 1.0])
    @test H == fill(2.0 + 0.0im, 1, 1, 2)
end

@testset "poles on the grid" begin
    integrator = transfer_function(zeros(1, 1), [1], [1], 0, s)
    H = frequency_response(integrator, s, [0.0, 1.0])
    @test !isfinite(H[1])
    @test H[2] ≈ -im
    oscillator = transfer_function([0 1; -1 0], [0, 1], [1, 0], 0, s)
    H = frequency_response(oscillator, s, [1.0, 2.0])
    @test !isfinite(H[1])
    @test H[2] ≈ -1 / 3
end

@testset "free parameters" begin
    G = transfer_function(fill(-1 / T, 1, 1), [1 / T], [1], 0, s)
    @test_throws ArgumentError frequency_response(G, s, [1.0])
    H = frequency_response(substitute(G, Dict(T => 2)), s, [0.5])
    @test H[1] ≈ 1 / (1 + im)
end

@testset "sizes" begin
    A, B, C, D = zeros(2, 2), zeros(2, 3), zeros(1, 2), zeros(1, 3)
    @test size(transfer_function(A, B, C, D, s)) == (1, 3)
    @test_throws DimensionMismatch transfer_function(zeros(2, 3), B, C, D, s)
    @test_throws DimensionMismatch transfer_function(A, zeros(3, 3), C, D, s)
    @test_throws DimensionMismatch transfer_function(A, zeros(3), C, D, s)
    @test_throws DimensionMismatch transfer_function(A, B, zeros(1, 3), D, s)
    @test_throws DimensionMismatch transfer_function(A, B, zeros(3), D, s)
    @test_throws DimensionMismatch transfer_function(A, B, C, zeros(3, 1), s)
    @test_throws DimensionMismatch transfer_function(A, B, C, zeros(3), s)
    @test_throws DimensionMismatch transfer_function(A, B, C, 0, s)
end
