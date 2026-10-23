# Tests for the public leadterm(expr, x, x0) (Task 4).

using Symbolics
using Symbolics: value, _nseries, series_coeff_at, series_valuation
using Test

@variables x a

# `leadterm` must be exported (not just Symbolics.leadterm) — use it unqualified below.

# compare a leadterm result (c, e): coeff matches (numerically or structurally) and exponent exact
function lead_is(expr, x, c_expected, e_expected; x0 = 0)
    c, e = leadterm(expr, x, x0)
    ce = isequal(value(value(c) - value(c_expected)), 0)
    return ce && e == Rational{Int}(e_expected)
end

@testset "leadterm" begin
    @testset "analytic" begin
        @test lead_is(sin(x), x, 1, 1)
        @test lead_is(cos(x), x, 1, 0)
        @test lead_is(exp(x), x, 1, 0)
        @test lead_is(x^3, x, 1, 3)
        @test lead_is(2x^2, x, 2, 2)
        @test lead_is(1 - cos(x), x, 1//2, 2)
        @test lead_is(x + x^2, x, 1, 1)
    end

    @testset "poles (negative exponent)" begin
        @test lead_is(exp(x)/x, x, 1, -1)
        @test lead_is(exp(x)/x^2, x, 1, -2)
        @test lead_is(1/x, x, 1, -1)
        @test lead_is(1/sin(x), x, 1, -1)
        @test lead_is(exp(x)/(x*(1+x)), x, 1, -1)
    end

    @testset "Puiseux (fractional exponent)" begin
        # x^(3//2): a pure fractional-power monomial with exponent > 1.
        # (leadterm(sqrt(x)) is covered via expand_series/_nseries in other testsets.)
        @test lead_is(x^(3//2), x, 1, 3//2)
    end

    @testset "Add-cancellation fallback" begin
        # sin(x) - x : x terms cancel, leading is -x^3/6
        @test lead_is(sin(x) - x, x, -1//6, 3)
        # cos(x) - 1 + x^2/2 = +x^4/24 - ...  (cos = 1 - x^2/2 + x^4/24)
        @test lead_is(cos(x) - 1 + x^2/2, x, 1//24, 4)
        # exp(x) - 1 - x : leading x^2/2
        @test lead_is(exp(x) - 1 - x, x, 1//2, 2)
    end

    @testset "shifted x0" begin
        # exp(x) about x0=1: leading term is exp(1) (SYMBOLIC, not 2.718) at exponent 0
        c, e = leadterm(exp(x), x, 1)
        @test e == 0
        @test !(c isa AbstractFloat)
        @test occursin("exp(1)", string(c))
        # (x-1) about x0=1: leading (x-1)^1
        @test lead_is(x - 1, x, 1, 1; x0 = 1)
    end

    @testset "symbolic / parametric coefficient" begin
        # a*x^2 : leading coeff a at exponent 2
        c, e = leadterm(a*x^2, x)
        @test e == 2
        @test isequal(value(value(c) - value(a)), 0)
    end

    @testset "complex rejection" begin
        # complex expressions have no real series in x -> clean ArgumentError
        @test_throws ArgumentError _nseries(exp(im*x), x, 3)
        @test_throws ArgumentError _nseries(im*x, x, 3)
        @test_throws ArgumentError _nseries(sin(im*x), x, 3)
    end

    @testset "essential singularity rejection" begin
        # A pole fed to exp/sin/cos/log is an essential singularity (or a log(x) term):
        # no Laurent/Puiseux series. Must reject cleanly (ArgumentError), NOT emit 1//0.
        @test_throws ArgumentError expand_series(exp(1/x), x, 3)
        @test_throws ArgumentError expand_series(exp(1/x^2), x, 3)
        @test_throws ArgumentError expand_series(sin(1/x), x, 3)
        @test_throws ArgumentError expand_series(cos(1/x), x, 3)
        @test_throws ArgumentError expand_series(log(x), x, 3)     # log(x) term at 0
        @test_throws ArgumentError expand_series(log(1/x), x, 3)   # -log(x) term
        # control: finite-nonzero inner (e0 == 0) is analytic and must still WORK
        s = expand_series(sin(x+1), x, 3)
        @test occursin("sin(1)", string(value(series_coeff_at(s, 0//1))))
        @test occursin("cos(1)", string(value(series_coeff_at(s, 1//1))))
    end

    @testset "Puiseux exact coefficients" begin
        # (2+x)^(1/2) = sqrt(2) + (1/4)sqrt(2) x - ...  (exact, irrational, not float).
        # NOTE: write the expected value as `Symbolics.term(sqrt, 2)`, NOT Julia's
        # `sqrt(2)` which eagerly evaluates to a Float64 and would defeat the check.
        s = _nseries((2+x)^(1//2), x, 2)
        c0 = value(series_coeff_at(s, 0//1))
        c1 = value(series_coeff_at(s, 1//1))
        sqrt2 = Symbolics.term(sqrt, 2)
        # coefficients must stay SYMBOLIC (sqrt(2)), not collapse to Float64
        @test !(c0 isa AbstractFloat)
        @test !(c1 isa AbstractFloat)
        # structural match against the symbolic expected forms
        @test isequal(Symbolics.unwrap(c0), Symbolics.unwrap(sqrt2))
        @test isequal(Symbolics.unwrap(c1), Symbolics.unwrap((1//4) * sqrt2))
        # perfect powers collapse exactly
        @test isequal(value(series_coeff_at(_nseries((4+x)^(1//2), x, 1), 0//1)), 2)
        # branch leading terms
        @test lead_is(1/sqrt(x), x, 1, -1//2)
        @test lead_is(x^(1//3) + x^(2//3), x, 1, 1//3)
    end

    @testset "expand_series public API" begin
        # exported (used unqualified here)
        # Taylor
        st = expand_series(sin(x), x, 6)
        @test isequal(value(series_coeff_at(st, 1//1)), 1)
        @test isequal(value(series_coeff_at(st, 3//1)), -1//6)
        @test isequal(value(series_coeff_at(st, 5//1)), 1//120)
        # Laurent (pole)
        sl = expand_series(exp(x)/x, x, 3)
        @test series_valuation(sl) == -1//1
        @test isequal(value(series_coeff_at(sl, -1//1)), 1)
        @test isequal(value(series_coeff_at(sl, 2//1)), 1//6)
        # Puiseux (branch point)
        sp = expand_series(sqrt(x), x, 2)
        @test series_valuation(sp) == 1//2
        @test isequal(value(series_coeff_at(sp, 1//2)), 1)
        # x0 shift keeps exp(1) symbolic
        ss = expand_series(exp(x), x, 2; x0 = 1)
        c = value(series_coeff_at(ss, 0//1))
        @test !(c isa AbstractFloat)
        @test occursin("exp(1)", string(c))
        # shifted pole: exp(x-1)/(x-1) about x0=1 is a simple pole
        ssp = expand_series(exp(x-1)/(x-1), x, 2; x0 = 1)
        @test series_valuation(ssp) == -1//1
        @test isequal(value(series_coeff_at(ssp, -1//1)), 1)
        # complex rejection propagates through the public API
        @test_throws ArgumentError expand_series(exp(im*x), x, 3)
    end

    @testset "pole + unit precision budget (Mathematica cross-check)" begin
        # 1/(x^2(1+x)) = x^-2 - x^-1 + 1 - x + x^2 - ...  To O(x^2) the x^1 = -1 term must
        # be present. Regression for the reciprocal maxk budget that undershot by `dv`
        # (the pole order) and dropped top terms. Mathematica: {1,-1,1,-1} from x^-2.
        s = expand_series(1/(x^2*(1+x)), x, 2)
        @test isequal(value(series_coeff_at(s, -2//1)), 1)
        @test isequal(value(series_coeff_at(s, -1//1)), -1)
        @test isequal(value(series_coeff_at(s, 0//1)), 1)
        @test isequal(value(series_coeff_at(s, 1//1)), -1)
        # exp(x)/x^2 to O(x^2) must carry through x^1 = 1/6
        s2 = expand_series(exp(x)/x^2, x, 2)
        @test isequal(value(series_coeff_at(s2, -2//1)), 1)
        @test isequal(value(series_coeff_at(s2, 1//1)), 1//6)
    end

    @testset "SymPy/Mathematica cross-validated" begin
        cvals(expr, ord, pairs) = begin
            s = expand_series(expr, x, ord)
            all(isequal(value(series_coeff_at(s, Rational{Int}(e))), c) for (e, c) in pairs)
        end
        # D1: 1/(x^3(1-2x))  = x^-3 + 2x^-2 + 4x^-1 + 8 + 16x   (pole x polynomial)
        @test cvals(1/(x^3*(1-2x)), 2, (-3=>1, -2=>2, -1=>4, 0=>8, 1=>16))
        # D2: 1/(x(1-x)^2)   = x^-1 + 2 + 3x + 4x^2
        @test cvals(1/(x*(1-x)^2), 3, (-1=>1, 0=>2, 1=>3, 2=>4))
        # D3: exp(sin(x))    = 1 + x + x^2/2 + 0 - x^4/8
        @test cvals(exp(sin(x)), 5, (0=>1, 1=>1, 2=>1//2, 3=>0, 4=>-1//8))
        # D4: log(cos(x))    = -x^2/2 - x^4/12 - x^6/45
        @test cvals(log(cos(x)), 7, (2=>-1//2, 4=>-1//12, 6=>-1//45))
        # D5: sin(x)/(1-cos(x)) = 2/x - x/6   (pole, transcendental denominator)
        @test cvals(sin(x)/(1-cos(x)), 3, (-1=>2, 1=>-1//6))
        # D6: tan(x)/x^2     = 1/x + x/3
        @test cvals(tan(x)/x^2, 3, (-1=>1, 1=>1//3))
        # D7: sqrt(sin(x))   = x^(1/2) - x^(5/2)/12   (Puiseux of a composite)
        @test cvals(sqrt(sin(x)), 3, (1//2=>1, 5//2=>-1//12))
        # D9: x^(2/3)(1+x)   = x^(2/3) + x^(5/3)
        @test cvals(x^(2//3)*(1+x), 3, (2//3=>1, 5//3=>1))
        # leadterm set (D13..D16)
        @test lead_is(tan(x)-sin(x), x, 1//2, 3)
        @test lead_is(1/x^2 - 1/x, x, 1, -2)
        @test lead_is(1-cos(x), x, 1//2, 2)
        @test lead_is(exp(x)-1-x-x^2/2, x, 1//6, 3)
    end
end
