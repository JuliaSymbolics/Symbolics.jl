# Unit tests for the TruncatedSeries core: hand-built series, no expression expansion.

using Symbolics
using Symbolics: TruncatedSeries, _series_monomial, _series_valuation, _ts_clean,
                 _series_coeff_at, _series_scale, SeriesDict
using Test

@variables a b

# build a series from `exponent => coefficient` pairs, truncated at `order`
_truncated_series(pairs, order) = TruncatedSeries(
    _ts_clean(Dict(Rational{Int}(e) => c for (e, c) in pairs), Rational{Int}(order)),
    Rational{Int}(order))

# tolerant zero check for symbolic cancellation
_ts_or_zero(c) = (c == 0) || isequal(Symbolics.value(c), 0)

@testset "TruncatedSeries core" begin
    @testset "construction + cleaning" begin
        s = _truncated_series([0 => 1, 1 => 2, 2 => 0, 3 => 4], 3)
        @test _series_coeff_at(s, 0) == 1
        @test _series_coeff_at(s, 1) == 2
        @test _series_coeff_at(s, 2) == 0   # zero coeff dropped
        @test _series_coeff_at(s, 3) == 0   # truncated at order 3
        @test s.order == 3
    end

    @testset "valuation" begin
        @test _series_valuation(_truncated_series([0 => 1, 1 => 2], 3)) == 0
        @test _series_valuation(_truncated_series([-2 => 1, 0 => 2], 3)) == -2   # Laurent
        @test _series_valuation(_truncated_series([1//2 => 1], 3)) == 1//2       # Puiseux
        @test _series_valuation(_truncated_series(SeriesDict(), 3)) == 3  # empty -> order
    end

    @testset "addition" begin
        s1 = _truncated_series([0 => 1, 1 => 2, 2 => 3], 4)
        s2 = _truncated_series([0 => 10, 1 => 20], 4)
        s = s1 + s2
        @test _series_coeff_at(s, 0) == 11
        @test _series_coeff_at(s, 1) == 22
        @test _series_coeff_at(s, 2) == 3
        @test s.order == 4

        s3 = _truncated_series([0 => 1, 1 => 5], 4)
        s4 = _truncated_series([0 => -1, 1 => 5], 4)
        sc = s3 + s4
        @test _series_coeff_at(sc, 0) == 0   # cancelled
        @test _series_coeff_at(sc, 1) == 10

        lo = _truncated_series([0 => 1, 1 => 1, 2 => 1], 2)
        hi = _truncated_series([0 => 1, 1 => 1, 2 => 1], 5)
        @test (lo + hi).order == 2   # order is the min
        @test _series_coeff_at(lo + hi, 2) == 0
    end

    @testset "multiplication" begin
        s = _truncated_series([0 => 1, 1 => 1], 3) * _truncated_series([0 => 1, 1 => 1], 3)
        @test _series_coeff_at(s, 0) == 1   # (1+x)^2 = 1 + 2x + x^2
        @test _series_coeff_at(s, 1) == 2
        @test _series_coeff_at(s, 2) == 1

        lau = _series_monomial(1, -1, 3)    # (x^-1)(1+x+x^2) = x^-1 + 1 + x
        tay = _truncated_series([0 => 1, 1 => 1, 2 => 1], 3)
        p = lau * tay
        @test _series_coeff_at(p, -1) == 1
        @test _series_coeff_at(p, 0) == 1
        @test _series_coeff_at(p, 1) == 1
        @test _series_valuation(p) == -1
        @test p.order == 2   # min(3+0, 3+(-1))

        ph = _series_monomial(1, 1//2, 2)   # x^(1/2) * x^(1/2) = x
        sq = ph * ph
        @test _series_coeff_at(sq, 1) == 1
        @test _series_valuation(sq) == 1
    end

    @testset "symbolic coefficients" begin
        s = _truncated_series([0 => a, 1 => 1], 3) * _truncated_series([0 => b, 1 => 1], 3)
        @test isequal(Symbolics.value(_series_coeff_at(s, 0)), Symbolics.value(a * b))
        @test isequal(Symbolics.value(_series_coeff_at(s, 1) - (a + b)), 0)
        @test _series_coeff_at(s, 2) == 1

        z = _truncated_series([0 => a], 2) + _truncated_series([0 => -a], 2)
        @test _ts_or_zero(_series_coeff_at(z, 0))   # a + (-a) drops
    end

    @testset "_series_scale" begin
        s = _truncated_series([0 => 1, 1 => 2, -1 => 3], 3)
        t = _series_scale(5, s)
        @test _series_coeff_at(t, -1) == 15
        @test _series_coeff_at(t, 0) == 5
        @test _series_coeff_at(t, 1) == 10
        @test _series_scale(0, s) |> isempty
    end
end
