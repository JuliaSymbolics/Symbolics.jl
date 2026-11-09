# Tests for _nseries: cross-checks vs taylor_coeff (analytic) and hand-known series.

using Symbolics
using Symbolics: _nseries, _series_coeff_at, _series_valuation, value
using Test

@variables x a b

# coeff of x^e from _nseries, as a plain value
nc(expr, e, order) = value(_series_coeff_at(_nseries(expr, x, order), Rational{Int}(e)))

# compare _nseries coeffs 0..hi against taylor_coeff
function matches_taylor(expr, hi)
    s = _nseries(expr, x, hi + 1)
    tc = Symbolics.taylor_coeff(expr, x, 0:hi; rationalize=true)
    all(isequal(value(_series_coeff_at(s, k)), value(tc[k+1])) for k in 0:hi)
end

@testset "_nseries series expansion" begin
    @testset "variable, const, polynomial" begin
        @test nc(x, 1, 3) == 1
        @test nc(x, 0, 3) == 0
        @test nc(3, 0, 3) == 3
        @test nc(1 + x + x^2, 2, 5) == 1
        @test nc(2x, 1, 3) == 2
    end

    @testset "exp / sin / cos vs taylor_coeff" begin
        @test matches_taylor(exp(x), 8)
        @test matches_taylor(sin(x), 8)
        @test matches_taylor(cos(x), 8)
        @test matches_taylor(exp(x) + sin(x), 7)
        @test matches_taylor(exp(x) * cos(x), 7)
        @test matches_taylor(sin(x^2), 8)       # vanishing inner arg of degree 2
        @test matches_taylor(exp(x^2), 8)
    end

    @testset "log vs taylor_coeff" begin
        @test matches_taylor(log(1 + x), 7)
        @test matches_taylor(log(1 - x), 7)
        @test matches_taylor(log(1 + x^2), 7)
    end

    @testset "integer powers / sqrt (Puiseux) vs taylor_coeff" begin
        @test matches_taylor((1 + x)^3, 5)
        @test matches_taylor((1 + x)^(-1), 6)     # geometric series, isdiv? ispow(-1)
        @test matches_taylor(sqrt(1 + x), 6)      # (1+x)^(1/2), analytic at 0
        # sqrt(x) itself is Puiseux: leading x^(1/2)
        @test nc(sqrt(x), 1//2, 3) == 1
        @test _series_valuation(_nseries(sqrt(x), x, 3)) == 1//2
    end

    @testset "Laurent poles (true poles, where taylor_coeff errors)" begin
        # exp(x)/x = x^-1 + 1 + x/2 + x^2/6 + ...
        @test nc(exp(x)/x, -1, 3) == 1
        @test nc(exp(x)/x, 0, 3) == 1
        @test nc(exp(x)/x, 1, 3) == 1//2
        @test nc(exp(x)/x, 2, 3) == 1//6
        @test _series_valuation(_nseries(exp(x)/x, x, 3)) == -1

        # exp(x)/x^2 = x^-2 + x^-1 + 1/2 + ...
        @test nc(exp(x)/x^2, -2, 2) == 1
        @test nc(exp(x)/x^2, -1, 2) == 1
        @test nc(exp(x)/x^2, 0, 2) == 1//2

        # 1/x
        @test nc(1/x, -1, 2) == 1
        @test nc(1/x, 0, 2) == 0

        # 1/sin(x) = x^-1 + x/6 + 7x^3/360 + ...
        @test nc(1/sin(x), -1, 4) == 1
        @test nc(1/sin(x), 1, 4) == 1//6
        @test nc(1/sin(x), 3, 4) == 7//360

        # non-monomial denominator: exp(x)/(x(1+x)) = x^-1 + 0 + x/2 - x^2/3 + ...
        @test nc(exp(x)/(x*(1+x)), -1, 3) == 1
        @test nc(exp(x)/(x*(1+x)), 0, 3) == 0
        @test nc(exp(x)/(x*(1+x)), 1, 3) == 1//2
        @test nc(exp(x)/(x*(1+x)), 2, 3) == -1//3
    end

    @testset "removable singularity" begin
        # sin(x)/x = 1 - x^2/6 + x^4/120
        @test nc(sin(x)/x, 0, 5) == 1
        @test nc(sin(x)/x, 2, 5) == -1//6
        @test nc(sin(x)/x, 4, 5) == 1//120
        # (1-cos(x))/x^2 = 1/2 - x^2/24 + ...
        @test nc((1-cos(x))/x^2, 0, 4) == 1//2
        @test nc((1-cos(x))/x^2, 2, 4) == -1//24
    end

    @testset "symbolic coefficients / free parameter" begin
        # a + sin(x)/x = (1+a) - x^2/6
        @test isequal(value(nc(a + sin(x)/x, 0, 3) - (1 + a)), 0)
        @test nc(a + sin(x)/x, 2, 3) == -1//6
        # a*exp(x)/x : pole with a parametric coefficient
        @test isequal(value(nc(a*exp(x)/x, -1, 2) - a), 0)
    end

    @testset "exact irrational coefficients (#1299)" begin
        # sin(x+1) coefficients (sin(1), cos(1), …) must stay exact-symbolic, never floats.
        s = _nseries(sin(x+1), x, 4)
        c0 = _series_coeff_at(s, 0//1)
        c1 = _series_coeff_at(s, 1//1)
        c2 = _series_coeff_at(s, 2//1)
        c3 = _series_coeff_at(s, 3//1)
        for c in (c0, c1, c2, c3)
            @test !(c isa AbstractFloat)
            @test c isa Symbolics.SymbolicT
        end
        @test occursin("sin(1)", string(c0))
        @test occursin("cos(1)", string(c1))
        @test occursin("sin(1)", string(c2)) && occursin("1//2", string(c2))
        @test occursin("cos(1)", string(c3)) && occursin("1//6", string(c3))
    end

    @testset "regression: Puiseux reciprocal (integer-grid bug)" begin
        # 1/(1+sqrt(x)) = 1 - x^(1/2) + x - x^(3/2) + …  (was silently 1)
        @test nc(1/(1+sqrt(x)), 0, 2) == 1
        @test nc(1/(1+sqrt(x)), 1//2, 2) == -1
        @test nc(1/(1+sqrt(x)), 1, 2) == 1
        @test nc(1/(1+sqrt(x)), 3//2, 2) == -1
        @test _series_valuation(_nseries(1/(1+sqrt(x)), x, 2)) == 0//1
        # cube-root grid
        cbrt_x = Symbolics.term(^, Symbolics.unwrap(x), 1//3)
        @test nc(1/(1+cbrt_x), 1//3, 1) == -1
        @test nc(1/(1+cbrt_x), 2//3, 1) == 1
        # precision-budget fix: 1/(x^2(1+x)) keeps its top terms
        @test nc(1/(x^2*(1+x)), -2, 2) == 1
        @test nc(1/(x^2*(1+x)), -1, 2) == -1
        @test nc(1/(x^2*(1+x)), 0, 2) == 1
        @test nc(1/(x^2*(1+x)), 1, 2) == -1
    end

    @testset "regression: adaptive leading-exponent probe (fixed-8 bug)" begin
        # 1/(sin(x^3)-x^3): den leads at x^9, so 1/den leads at x^-9 (was 'bad valuation 8//1')
        @test nc(1/(sin(x^3)-x^3), -9, -7) == -6
        @test _series_valuation(_nseries(1/(sin(x^3)-x^3), x, -7)) == -9//1
    end

    @testset "regression: trig of non-analytic finite argument (1//0 leak)" begin
        # sin(1+sqrt(x)) via angle addition (was (1//0)*x); coeffs stay exact-symbolic
        s = _nseries(sin(1+sqrt(x)), x, 2)
        c_h = _series_coeff_at(s, 1//2)   # cos(1)
        c_1 = _series_coeff_at(s, 1//1)   # -sin(1)/2
        @test !(value(c_h) isa AbstractFloat)
        @test occursin("cos(1)", string(c_h))
        @test occursin("sin(1)", string(c_1)) && occursin("1//2", string(c_1))
        @test !any(occursin("1//0", string(value(c))) for c in values(s.coeffs))
    end

    @testset "regression: parameter cancellation leaves no phantom zero (Bug 3 guard)" begin
        # x^1 coefficient must collapse to a true zero; leading exponent is 2, not a phantom
        e = (a + b)*x - a*x - b*x + x^2
        @test _series_valuation(_nseries(e, x, 3)) == 2//1
        @test nc(e, 2, 3) == 1
        @test value(_series_coeff_at(_nseries(e, x, 3), 1//1)) == 0
    end

    @testset "regression: high-order factorial overflow (big coefficients)" begin
        # Int64 coeffs overflowed ~order 13 in composed functions; coeffs now widen to BigInt
        # as needed and narrow back when they fit. Expected values are self-evident 1/n!.
        @test nc(exp(x), 20, 21) == 1 // factorial(big(20))
        @test nc(exp(x), 25, 26) == 1 // factorial(big(25))
        s = _nseries(exp(sqrt(x)), x, 30)   # Puiseux + high order: no overflow
        @test _series_valuation(s) == 0//1
        @test nc(exp(sqrt(x)), 10, 30) == 1 // factorial(big(20))
        @test _series_coeff_at(_nseries(exp(x), x, 5), 2//1) isa Rational{Int}   # narrows back
    end

    @testset "regression: unary functions of non-analytic arguments (1//0 leak)" begin
        # Reference values from SymPy series(); these previously returned (1//0)*x.
        @test nc(tan(sqrt(x)), 1//2, 2) == 1 && nc(tan(sqrt(x)), 3//2, 2) == 1//3
        @test nc(sinh(sqrt(x)), 1//2, 2) == 1 && nc(sinh(sqrt(x)), 3//2, 2) == 1//6
        s = _nseries(atan(1 + sqrt(x)), x, 2)   # atan(1) + x^(1/2)/2 - x/4 + x^(3/2)/12
        @test isequal(value(_series_coeff_at(s, 0)), value(atan(Num(1))))
        @test [value(_series_coeff_at(s, e)) for e in (1//2, 1, 3//2)] == [1//2, -1//4, 1//12]
        @test nc(tanh(x^2), 6, 7) == -1//3
        @test nc(asin(x/2), 3, 4) == 1//48
        @test matches_taylor(tan(x), 6) && matches_taylor(atan(x), 6)
        @test isequal(nc(tan(a*x), 3, 4), value((1//3) * a^3))
        @test nc(tan(x)/x^2, -1, 3) == 1 && nc(tan(x)/x^2, 1, 3) == 1//3
        @test nc(atan(x, 1 + x), 2, 3) == -1                    # two-argument term
    end

    @testset "regression: non-analytic unary functions are rejected" begin
        @test_throws ArgumentError _nseries(asin(1 + x), x, 2)   # branch point of asin at 1
        @test_throws ArgumentError _nseries(atan(1/x), x, 2)     # pole inside the function
        @test_throws ArgumentError _nseries(abs(x), x, 2)        # kink: was silently x
        @test_throws ArgumentError _nseries(sign(x), x, 2)       # jump: was silently 0
        @test nc(abs(1 + x), 1, 3) == 1                          # smooth away from the kink
    end
end
