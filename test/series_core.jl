# Unit tests for the TruncatedSeries core (Task 2). Hand-built series only; no
# expression expansion yet (that is _nseries, a later task).

using Symbolics
using Symbolics: TruncatedSeries, truncated_series, series_monomial, series_valuation,
                 series_coeff_at, series_scale
using Test

@variables a b

# tolerant zero check for symbolic cancellation result
_ts_or_zero(c) = (c == 0) || isequal(Symbolics.value(c), 0)

@testset "TruncatedSeries core" begin
    @testset "construction + cleaning" begin
        s = truncated_series([0 => 1, 1 => 2, 2 => 0, 3 => 4], 3)
        # order 3 drops the x^3 term; zero coeff at x^2 dropped
        @test series_coeff_at(s, 0) == 1
        @test series_coeff_at(s, 1) == 2
        @test series_coeff_at(s, 2) == 0   # absent
        @test series_coeff_at(s, 3) == 0   # truncated
        @test s.order == 3
    end

    @testset "valuation" begin
        @test series_valuation(truncated_series([0 => 1, 1 => 2], 3)) == 0
        @test series_valuation(truncated_series([-2 => 1, 0 => 2], 3)) == -2   # Laurent
        @test series_valuation(truncated_series([1//2 => 1], 3)) == 1//2       # Puiseux
        @test series_valuation(truncated_series(Dict{Rational{Int},Any}(), 3)) == 3  # empty -> order
    end

    @testset "addition" begin
        s1 = truncated_series([0 => 1, 1 => 2, 2 => 3], 4)
        s2 = truncated_series([0 => 10, 1 => 20], 4)
        s = s1 + s2
        @test series_coeff_at(s, 0) == 11
        @test series_coeff_at(s, 1) == 22
        @test series_coeff_at(s, 2) == 3
        @test s.order == 4

        # cancellation drops the term
        s3 = truncated_series([0 => 1, 1 => 5], 4)
        s4 = truncated_series([0 => -1, 1 => 5], 4)
        sc = s3 + s4
        @test series_coeff_at(sc, 0) == 0   # cancelled, absent
        @test series_coeff_at(sc, 1) == 10

        # order is the min of the two
        lo = truncated_series([0 => 1, 1 => 1, 2 => 1], 2)
        hi = truncated_series([0 => 1, 1 => 1, 2 => 1], 5)
        @test (lo + hi).order == 2
        @test series_coeff_at(lo + hi, 2) == 0  # beyond order 2
    end

    @testset "multiplication" begin
        # (1 + x)(1 + x) = 1 + 2x + x^2
        s = truncated_series([0 => 1, 1 => 1], 3) * truncated_series([0 => 1, 1 => 1], 3)
        @test series_coeff_at(s, 0) == 1
        @test series_coeff_at(s, 1) == 2
        @test series_coeff_at(s, 2) == 1

        # Laurent * Taylor: (x^-1)(1 + x + x^2) = x^-1 + 1 + x
        lau = series_monomial(1, -1, 3)
        tay = truncated_series([0 => 1, 1 => 1, 2 => 1], 3)
        p = lau * tay
        @test series_coeff_at(p, -1) == 1
        @test series_coeff_at(p, 0) == 1
        @test series_coeff_at(p, 1) == 1
        @test series_valuation(p) == -1

        # precision bookkeeping: val(lau) = -1, so order = min(3 + 0, 3 + (-1)) = 2
        @test p.order == 2

        # Puiseux * Puiseux: x^(1/2) * x^(1/2) = x^1
        ph = series_monomial(1, 1//2, 2)
        sq = ph * ph
        @test series_coeff_at(sq, 1) == 1
        @test series_valuation(sq) == 1
    end

    @testset "symbolic coefficients" begin
        # (a + x)(b + x) = a*b + (a+b)x + x^2, coefficients are symbolic
        s = truncated_series([0 => a, 1 => 1], 3) * truncated_series([0 => b, 1 => 1], 3)
        @test isequal(Symbolics.value(series_coeff_at(s, 0)), Symbolics.value(a * b))
        @test isequal(Symbolics.value(series_coeff_at(s, 1) - (a + b)), 0)
        @test series_coeff_at(s, 2) == 1

        # symbolic cancellation: (a)x^0 + (-a)x^0 should drop (structural)
        z = truncated_series([0 => a], 2) + truncated_series([0 => -a], 2)
        # -a + a is structurally 0 after value() in most cases; accept either absent or zero
        @test _ts_or_zero(series_coeff_at(z, 0))
    end

    @testset "series_scale" begin
        s = truncated_series([0 => 1, 1 => 2, -1 => 3], 3)
        t = series_scale(5, s)
        @test series_coeff_at(t, -1) == 15
        @test series_coeff_at(t, 0) == 5
        @test series_coeff_at(t, 1) == 10
        @test series_scale(0, s) |> isempty
    end
end
