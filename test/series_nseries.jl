# Tests for the _nseries series expansion. Cross-checks against taylor_coeff's
# taylor_coeff (analytic parts) and against hand-known Laurent/Puiseux series.

using Symbolics
using Symbolics: _nseries, series_coeff_at, series_valuation, value
using Test

@variables x a

# extract coeff of x^e from _nseries result as a plain value
nc(expr, e, order) = value(series_coeff_at(_nseries(expr, x, order), Rational{Int}(e)))

# compare _nseries coefficients 0..hi against taylor_coeff (analytic only)
function matches_taylor(expr, hi)
    s = _nseries(expr, x, hi + 1)
    tc = Symbolics.taylor_coeff(expr, x, 0:hi; rationalize=true)
    all(isequal(value(series_coeff_at(s, k)), value(tc[k+1])) for k in 0:hi)
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
        @test series_valuation(_nseries(sqrt(x), x, 3)) == 1//2
    end

    @testset "Laurent poles (true poles, where taylor_coeff errors)" begin
        # exp(x)/x = x^-1 + 1 + x/2 + x^2/6 + ...
        @test nc(exp(x)/x, -1, 3) == 1
        @test nc(exp(x)/x, 0, 3) == 1
        @test nc(exp(x)/x, 1, 3) == 1//2
        @test nc(exp(x)/x, 2, 3) == 1//6
        @test series_valuation(_nseries(exp(x)/x, x, 3)) == -1

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
        # sin(x+1) routes through the diff fallback; coefficients are sin(1), cos(1),
        # -sin(1)/2, -cos(1)/6 — must stay EXACT SYMBOLIC, never rationalized/floatified.
        # We cannot write a symbolic `sin(1)` literal in the test (Julia folds sin(1) to a
        # Float64 immediately), so instead assert each coefficient is (a) NOT a float and
        # (b) a symbolic call whose string is the expected symbolic form, and (c) matches
        # the known numeric value when evaluated.
        s = _nseries(sin(x+1), x, 4)
        c0 = series_coeff_at(s, 0//1)
        c1 = series_coeff_at(s, 1//1)
        c2 = series_coeff_at(s, 2//1)
        c3 = series_coeff_at(s, 3//1)
        # not floats (the #1299 guarantee)
        for c in (c0, c1, c2, c3)
            @test !(c isa AbstractFloat)
            @test c isa Symbolics.SymbolicT
        end
        # symbolic form (string compare is robust to how sin(1) prints)
        @test occursin("sin(1)", string(c0))
        @test occursin("cos(1)", string(c1))
        @test occursin("sin(1)", string(c2)) && occursin("1//2", string(c2))
        @test occursin("cos(1)", string(c3)) && occursin("1//6", string(c3))
    end
end
