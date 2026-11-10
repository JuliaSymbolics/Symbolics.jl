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
        # reciprocal recurrence and binomial coefficients overflowed Int64 (SymPy values)
        @test nc(x/sin(x), 20, 22) == 91546277357 // big(401428831349145600000)
        @test nc((1 + x)^(1//2), 39, 40) ==
              big(22091732848375837675) // big(18889465931478580854784)
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

    @testset "regression: integer powers keep the requested precision" begin
        s = _nseries((1/x + 1)^3, x, 3)        # was truncated at O(x^1)
        @test s.order == 3//1 && nc((1/x + 1)^3, 0, 3) == 1 && nc((1/x + 1)^3, -3, 3) == 1
        s = _nseries((x + x^2)^3, x, 5)        # returned terms through x^6 + O(x^7)
        @test s.order == 5//1 && maximum(keys(s.coeffs)) == 4//1
        @test nc((1/sin(x))^2, 2, 3) == 1//15  # SymPy: x^-2 + 1/3 + x^2/15
        @test nc(sin(x)^2, 6, 7) == 2//45
        @test _nseries(x^6, x, 5).order == 5//1 && isempty(_nseries(x^6, x, 5))
        @test nc((1 + x)^1000, 9, 10) == binomial(big(1000), 9)
    end

    @testset "regression: factors and numerators with poles keep the requested precision" begin
        s = _nseries((1/x^2 + x)/(1 + x), x, 2)    # SymPy: x^-2 - 1/x + 1 + O(x^2)
        @test s.order == 2//1 && nc((1/x^2 + x)/(1 + x), 0, 2) == 1
        s = _nseries(exp(x)*x^-2, x, 2)            # SymPy: x^-2 + 1/x + 1/2 + x/6 + O(x^2)
        @test s.order == 2//1 && nc(exp(x)*x^-2, 1, 2) == 1//6
        @test nc(sin(x)*x^-3, 0, 2) == -1//6
        @test _nseries(x^-2*sqrt(sin(x)), x, 0).order == 0//1
        s = _nseries(x^5*sqrt(sin(x)), x, 3)       # SymPy: O(x^3)
        @test s.order == 3//1 && isempty(s)
        @test _nseries((sin(x)^2 + cos(x)^2 - 1)/(1 + x), x, 3).order == 3//1
    end

    @testset "regression: fractional power truncated below its leading term" begin
        for (f, o) in ((sqrt(sin(x)), 0//1), (sqrt(sin(x)), 1//2), (x^(1//3), 1//3), (sin(x)^(-1//2), -1//1))
            s = _nseries(f, x, o)                  # threw "base leading coefficient zero"
            @test s.order == o && isempty(s)
        end
        @test string(expand_series(sqrt(sin(x)) + 1, x, 0)) == "O(x^0)"
    end

    @testset "regression: exact constant roots (Float64 check gave wrong values)" begin
        @test !(nc(sqrt(big(10)^18 + 1 + x), 0, 2) isa Number)      # was 1000000000
        @test nc(sqrt(big(10)^18 + x), 0, 2) == 10^9
        @test nc(sqrt(1//4 + x), 0, 2) == 1//2 && nc(sqrt(1//4 + x), 1, 2) == 1
        @test nc((8//27 + x)^(1//3), 1, 2) == 3//4
        @test nc((big(10)^10 + x)^(-3), 1, 2) == -3 // big(10)^40     # Int64 overflow
        @test Symbolics._pow_const(big(3)^90, 1//3) == big(3)^30
        @test !(Symbolics._pow_const(-8, 1//3) isa Number)       # principal root is complex
    end

    @testset "pole rewrites and real roots at a zero of the argument" begin
        # SymPy: cot = 1/x - x/3 - x^3/45, csc = 1/x + x/6 + 7x^3/360, coth/csch alternate
        for (f, c) in ((cot(x), (1, -1//3, -1//45)), (csc(x), (1, 1//6, 7//360)),
                       (coth(x), (1, 1//3, -1//45)), (csch(x), (1, -1//6, 7//360)))
            @test (nc(f, -1, 4), nc(f, 1, 4), nc(f, 3, 4)) == c
        end
        @test nc(cot(2x), -1, 2) == 1//2 && nc(cot(2x), 1, 2) == -2//3
        @test nc(csc(sin(x)), 1, 2) == 1//3
        @test nc(cot(x)^2, 0, 1) == -2//3
        @test string(nc(cot(1 + x), 0, 2)) == "cot(1)"         # unchanged at a nonzero point
        @test string(expand_series(cbrt(-x + x^2), x, 2)) == "-x^(1//3) + (1//3)*x^(4//3) + O(x^2)"
        @test nc(cbrt(sin(x)), 7//3, 3) == -1//18
        @test nc(cbrt(1/x), -1//3, 2) == 1
        @test nc(cbrt(-8 + x), 0, 2) == -2
        @test_throws ArgumentError expand_series(cot(1/x), x, 2)
        @test nc(cot(0.5x), -1, 2) == 2.0 && nc(cot(0.5x), 1, 2) ≈ -1/6
    end

    @testset "Float inputs (SymPy / Mathematica values)" begin
        isf(c) = value(c) isa Float64
        @test nc(exp(0.1x), 1, 3) ≈ 0.1 && nc(exp(0.1x), 2, 3) ≈ 0.005 && isf(nc(exp(0.1x), 2, 3))
        @test nc(sin(x) + 0.5x, 1, 4) == 1.5 && nc(sin(x) + 0.5x, 3, 4) == -1//6
        @test nc(sqrt(2.0 + x), 0, 2) ≈ sqrt(2.0) && nc(sqrt(2.0 + x), 1, 2) ≈ 0.3535533905932738
        @test nc(log(1.5 + x), 0, 2) ≈ log(1.5) && nc(1/(2.0 + x), 2, 3) == 0.125
        @test nc(sin(0.5 + x), 2, 3) ≈ -0.2397127693021015
        @test nc(log(cos(0.5x)), 4, 5) ≈ -0.005208333333333334
        @test nc(cos(0.5x)/sin(0.5x), -1, 2) == 2.0 && nc(cos(0.5x)/sin(0.5x), 1, 2) ≈ -1/6
        @test nc(1/(2.0 + sqrt(x)), 1//2, 2) == -0.25
        @test nc(1 + 0.5/x, -1, 3) == 0.5
        @test nc((1 + x)^0.5, 2, 3) == -0.125                  # float exponent, rationalized
        @test nc(x^1.5, 3//2, 3) == 1 && nc(x^0.1, 1//10, 1) == 1 && nc(sqrt(x)^3.0, 3//2, 3) == 1
        @test nc((-2.0 + x)^3.0, 1, 2) == 12
        @test string(expand_series(exp(x), x, 2; x0 = 0.5)) ==
              "1.6487212707001282 + 1.6487212707001282*(x - 0.5) + O((x - 0.5)^2)"
        for (f, o) in ((exp(0.1x), 30), (0.5a*exp(x), 25), (exp(0.5a*x), 24), (x/sin(0.5x), 24))
            s = _nseries(f, x, o)                              # no BigFloat promotion
            @test !any(c -> occursin("BigFloat", string(typeof(value(c)))) ||
                            occursin(r"\d{25,}", string(c)), values(s.coeffs))
        end
        @test nc(Float32(0.5)*x + 1, 1, 2) === 0.5 && nc(big"0.5"*x, 1, 2) === 0.5
    end

    @testset "printing: variable, expansion point, fractional exponents" begin
        @variables t
        @test string(expand_series(sin(t), t, 4)) == "t + (-1//6)*t^3 + O(t^4)"
        @test string(expand_series(1/x, x, 3; x0 = 2)) ==
              "1//2 + (-1//4)*(x - 2) + (1//8)*(x - 2)^2 + O((x - 2)^3)"
        @test string(expand_series(exp(x), x, 2; x0 = -1)) ==
              "exp(-1) + (exp(-1))*(x + 1) + O((x + 1)^2)"
        @test string(expand_series(1/(1 + sqrt(x)), x, 2)) ==
              "1 - x^(1//2) + x - x^(3//2) + O(x^2)"
        @test string(expand_series(exp(x)/x, x, 2)) == "x^-1 + 1 + (1//2)*x + O(x^2)"
        @test string(expand_series(sin(x), x, 0)) == "O(x^0)"
    end
end
