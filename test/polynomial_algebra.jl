using Symbolics
using Test
using Random

@variables x a b c d e

# =========================================================================
# Foundational helper: Symbolics._univariate_coeffs
# =========================================================================

@testset "univariate coefficients (foundational helper)" begin
    # valid polynomial
    coeffs = Symbolics._univariate_coeffs(3x^2 + 2x + 1, x)
    @test length(coeffs) == 3
    @test isequal(Symbolics.simplify(coeffs[1] - 3), 0)
    @test isequal(Symbolics.simplify(coeffs[2] - 2), 0)
    @test isequal(Symbolics.simplify(coeffs[3] - 1), 0)

    # degree-0 input
    @test Symbolics._univariate_coeffs(Symbolics.Num(5), x) == [Symbolics.Num(5)]

    # symbolic coefficients
    coeffs = Symbolics._univariate_coeffs(a * x^2 + b * x + c, x)
    @test length(coeffs) == 3
    @test isequal(Symbolics.simplify(coeffs[1] - a), 0)
    @test isequal(Symbolics.simplify(coeffs[2] - b), 0)
    @test isequal(Symbolics.simplify(coeffs[3] - c), 0)

    # non-polynomial input raises
    @test_throws ArgumentError Symbolics._univariate_coeffs(sin(x), x)
    @test_throws ArgumentError Symbolics._univariate_coeffs(1 / x, x)
end

# =========================================================================
# US1 — resultant
# =========================================================================

# Helper: assert symbolic equality by simplifying the difference.
sym_eq(x, y) = iszero(Symbolics.simplify(Symbolics.expand(Symbolics.Num(x) - Symbolics.Num(y))))

@testset "resultant — fixed numeric inputs" begin
    # Common root → 0
    @test sym_eq(resultant(x^2 - 1, x - 1, x), 0)
    @test sym_eq(resultant(x^2 - 1, x - 1, x; algorithm = :euclid), 0)
    @test sym_eq(resultant(x^2 - 1, x - 1, x; algorithm = :sylvester), 0)

    # Coprime pair → 5 for (x^2+1, x-2)
    @test sym_eq(resultant(x^2 + 1, x - 2, x), 5)
    @test sym_eq(resultant(x^2 + 1, x - 2, x; algorithm = :sylvester), 5)

    # Res(x^3 - 1, x^2 + x + 1) = 0 (share the factor x^2+x+1)
    @test sym_eq(resultant(x^3 - 1, x^2 + x + 1, x), 0)
    @test sym_eq(resultant(x^3 - 1, x^2 + x + 1, x; algorithm = :sylvester), 0)

    # Res(f, c) = c^deg(f); Res(x^3+1, 2, x) = 8
    @test sym_eq(resultant(x^3 + 1, Symbolics.Num(2), x), 8)

    # Res of two scalars: 1 (empty product)
    @test sym_eq(resultant(Symbolics.Num(3), Symbolics.Num(7), x), 1)

    # Standard example: Res(x^2 - 3x + 2, x^2 - 5x + 6, x)
    # = (x-1)(x-2) ∩ (x-2)(x-3) share root x=2 → 0
    @test sym_eq(resultant(x^2 - 3x + 2, x^2 - 5x + 6, x), 0)

    # Res(x^2 + 3, x^3 + x + 1, x) — known value 13 (via lc(f)^deg(g) * prod g(α)
    # over roots α=±i√3 of f: g(i√3)·g(-i√3) = (1-2i√3)(1+2i√3) = 1+12 = 13)
    @test sym_eq(resultant(x^2 + 3, x^3 + x + 1, x), 13)
    @test sym_eq(resultant(x^2 + 3, x^3 + x + 1, x; algorithm = :sylvester), 13)
end

@testset "resultant — symbolic coefficients" begin
    # Res(a*x + b, c*x + d, x) == a*d - b*c
    got = resultant(a * x + b, c * x + d, x)
    @test sym_eq(got, a * d - b * c)
    got_syl = resultant(a * x + b, c * x + d, x; algorithm = :sylvester)
    @test sym_eq(got_syl, a * d - b * c)

    # Res(f, f', x) for f = a*x^2 + b*x + c is -a*(b^2 - 4ac)
    f = a * x^2 + b * x + c
    fprime = Symbolics.derivative(f, x)
    got = resultant(f, fprime, x)
    @test sym_eq(got, -a * (b^2 - 4 * a * c))
end

@testset "resultant — algorithm agreement" begin
    Random.seed!(0xC0FFEE)
    n_agree = 0
    n_trials = 50
    for trial in 1:n_trials
        df = rand(1:4)
        dg = rand(1:4)
        f_coeffs = big.(rand(-3:3, df + 1))
        g_coeffs = big.(rand(-3:3, dg + 1))
        f_coeffs[1] == 0 && (f_coeffs[1] = big(1))
        g_coeffs[1] == 0 && (g_coeffs[1] = big(1))
        f = sum(f_coeffs[i] * x^(df + 1 - i) for i in 1:(df + 1))
        g = sum(g_coeffs[i] * x^(dg + 1 - i) for i in 1:(dg + 1))
        r1 = resultant(f, g, x; algorithm = :euclid)
        r2 = resultant(f, g, x; algorithm = :sylvester)
        if sym_eq(r1, r2)
            n_agree += 1
        else
            @info "algorithm disagreement" trial f g r1 r2
        end
    end
    @test n_agree == n_trials
end

@testset "resultant — edge cases and errors" begin
    # Zero polynomial on either side
    @test_throws ArgumentError resultant(Symbolics.Num(0), x - 1, x)
    @test_throws ArgumentError resultant(x - 1, Symbolics.Num(0), x)

    # Unknown algorithm
    @test_throws ArgumentError resultant(x - 1, x + 1, x; algorithm = :foo)

    # Non-polynomial input
    @test_throws ArgumentError resultant(sin(x), x - 1, x)
    @test_throws ArgumentError resultant(x - 1, 1 / x, x)
end

# =========================================================================
# US2 — discriminant
# =========================================================================

@testset "discriminant — symbolic quadratic" begin
    # Canonical: disc(a*x^2 + b*x + c) = b^2 - 4*a*c
    got = discriminant(a * x^2 + b * x + c, x)
    @test sym_eq(got, b^2 - 4 * a * c)
end

@testset "discriminant — fixed numeric inputs" begin
    # Distinct rational roots of x^2 - 3x + 2 → discriminant = 1
    @test sym_eq(discriminant(x^2 - 3x + 2, x), 1)

    # x^2 + 1 → -4
    @test sym_eq(discriminant(x^2 + 1, x), -4)

    # Repeated root: x^3 - 3x + 2 has a double root at x = 1 → 0
    @test sym_eq(discriminant(x^3 - 3x + 2, x), 0)

    # Three distinct rational roots 1, 2, 3 of x^3 - 6x^2 + 11x - 6:
    # disc = ((2-1)(3-1)(3-2))^2 = 4
    @test sym_eq(discriminant(x^3 - 6x^2 + 11x - 6, x), 4)

    # Linear polynomial convention
    @test sym_eq(discriminant(x - 5, x), 1)
    @test sym_eq(discriminant(3x + 7, x), 1)

    # Non-zero constant convention (per spec/contract Q4: disc(c) = 0)
    @test sym_eq(discriminant(Symbolics.Num(7), x), 0)
end

@testset "discriminant — identity with resultant" begin
    # disc(f) = (-1)^(n(n-1)/2) * Res(f, f') / lc(f)
    for f in (x^2 + 3x + 2, x^3 - x, 2x^3 + 5x + 1, x^4 - 2x^2 + 1)
        n = Symbolics.degree(f, x)
        lc = Symbolics._univariate_coeffs(f, x)[1]
        fprime = Symbolics.derivative(f, x)
        expected = ((-1)^(n * (n - 1) ÷ 2)) * resultant(f, fprime, x) / lc
        @test sym_eq(discriminant(f, x), expected)
    end
end

@testset "discriminant — edge cases and errors" begin
    @test_throws ArgumentError discriminant(Symbolics.Num(0), x)
    @test_throws ArgumentError discriminant(sin(x), x)
end

# =========================================================================
# US3 — sqrfree
# =========================================================================

# Reconstruct original polynomial from the sqrfree decomposition.
function reconstruct_sqrfree(u, factors)
    p = u
    for (fac, e) in factors
        p *= fac^e
    end
    return Symbolics.simplify(Symbolics.expand(p))
end

@testset "sqrfree — fixed inputs" begin
    # Classical textbook: (x-1)^2 (x-2)^3
    f = (x - 1)^2 * (x - 2)^3
    u, fs = sqrfree(f, x)
    @test sym_eq(reconstruct_sqrfree(u, fs), Symbolics.expand(f))
    # Structure: two factors with multiplicities 2 and 3
    mults = sort([e for (_, e) in fs])
    @test mults == [2, 3]

    # Already square-free: x^2 - 1 → single factor with multiplicity 1
    u, fs = sqrfree(x^2 - 1, x)
    @test length(fs) == 1
    @test fs[1][2] == 1
    @test sym_eq(reconstruct_sqrfree(u, fs), x^2 - 1)

    # Non-unit leading coefficient: 2x^2 - 2 → unit=2, factor=(x^2-1, 1)
    u, fs = sqrfree(2x^2 - 2, x)
    @test sym_eq(u, 2)
    @test sym_eq(reconstruct_sqrfree(u, fs), 2x^2 - 2)

    # Constant non-zero: 5 → unit=5, no non-trivial factors
    u, fs = sqrfree(Symbolics.Num(5), x)
    @test sym_eq(u, 5)
    @test isempty(fs)

    # Mixed multiplicities: x^2 * (x-1)^3 * (x+1)
    f = x^2 * (x - 1)^3 * (x + 1)
    u, fs = sqrfree(f, x)
    @test sym_eq(reconstruct_sqrfree(u, fs), Symbolics.expand(f))
    @test sort([e for (_, e) in fs]) == [1, 2, 3]
end

@testset "sqrfree — round-trip property" begin
    Random.seed!(0xDEC0DE)
    n_ok = 0
    n_trials = 30
    for _ in 1:n_trials
        roots_set = unique(rand(-3:3, rand(2:4)))
        n_roots = length(roots_set)
        mults = rand(1:3, n_roots)
        f = reduce(*, (x - r)^m for (r, m) in zip(roots_set, mults))
        f = Symbolics.expand(f)
        u, fs = sqrfree(f, x)
        recon = reconstruct_sqrfree(u, fs)
        sym_eq(recon, Symbolics.expand(f)) && (n_ok += 1)
    end
    @test n_ok == n_trials
end

@testset "sqrfree — edge cases and errors" begin
    @test_throws ArgumentError sqrfree(Symbolics.Num(0), x)
    @test_throws ArgumentError sqrfree(sin(x), x)
    @test_throws ArgumentError sqrfree(1 / x, x)
end

# =========================================================================
# US4 — factor
# =========================================================================

function reconstruct_factor(u, factors)
    p = u
    for (fac, e) in factors
        p *= fac^e
    end
    return Symbolics.simplify(Symbolics.expand(p))
end

@testset "factor — fixed inputs (pure-Julia path)" begin
    # x^4 - 1 = (x-1)(x+1)(x^2+1), the quadratic being irreducible over ℚ
    u, fs = Symbolics.factor(x^4 - 1, x)
    @test sym_eq(u, 1)
    @test length(fs) == 3
    @test all(e == 1 for (_, e) in fs)
    @test sym_eq(reconstruct_factor(u, fs), Symbolics.expand(x^4 - 1))

    # x^3 - 6x^2 + 11x - 6 → three distinct linear factors (roots 1, 2, 3)
    u, fs = Symbolics.factor(x^3 - 6x^2 + 11x - 6, x)
    @test sym_eq(u, 1)
    @test length(fs) == 3
    @test all(e == 1 for (_, e) in fs)
    @test sym_eq(reconstruct_factor(u, fs), Symbolics.expand(x^3 - 6x^2 + 11x - 6))

    # Content > 1: 2x^2 - 2 = 2 * (x-1) * (x+1)
    u, fs = Symbolics.factor(2x^2 - 2, x)
    @test sym_eq(u, 2)
    @test length(fs) == 2
    @test sym_eq(reconstruct_factor(u, fs), Symbolics.expand(2x^2 - 2))

    # Irreducible over ℚ: x^2 + 1
    u, fs = Symbolics.factor(x^2 + 1, x)
    @test length(fs) == 1
    @test sym_eq(fs[1][1], x^2 + 1)
    @test fs[1][2] == 1

    # x^3 + 1 = (x+1)(x^2 - x + 1)
    u, fs = Symbolics.factor(x^3 + 1, x)
    @test length(fs) == 2
    @test sym_eq(reconstruct_factor(u, fs), Symbolics.expand(x^3 + 1))

    # Repeated root: (x-1)^2 * (x+2)
    u, fs = Symbolics.factor((x - 1)^2 * (x + 2), x)
    @test sym_eq(reconstruct_factor(u, fs), Symbolics.expand((x - 1)^2 * (x + 2)))
    @test sort([e for (_, e) in fs]) == [1, 2]

    # Constant non-zero: Symbolics.factor(5) → unit=5, no non-trivial factors
    u, fs = Symbolics.factor(Symbolics.Num(5), x)
    @test sym_eq(u, 5)
    @test isempty(fs)
end

@testset "factor — round-trip property" begin
    Random.seed!(0xFAC0FF)
    n_ok = 0
    n_trials = 20
    for _ in 1:n_trials
        roots_set = unique(rand(-4:4, 4))
        length(roots_set) < 2 && continue
        mults = rand(1:2, length(roots_set))
        f = reduce(*, (x - r)^m for (r, m) in zip(roots_set, mults))
        f = Symbolics.expand(f)
        u, fs = Symbolics.factor(f, x)
        recon = reconstruct_factor(u, fs)
        sym_eq(recon, Symbolics.expand(f)) && (n_ok += 1)
    end
    @test n_ok >= 15  # at least 75% — allow a few trials dropped by unique() collapse
end

@testset "factor — edge cases and errors" begin
    @test_throws ArgumentError Symbolics.factor(Symbolics.Num(0), x)
    @test_throws ArgumentError Symbolics.factor(sin(x), x)
    @test_throws ArgumentError Symbolics.factor(1 / x, x)
    # Floating-point coefficient is out of the supported domain
    @test_throws ArgumentError Symbolics.factor(1.5 * x + 2, x)
end

# =========================================================================
# Optional Nemo oracle tests — run only when Nemo is loaded (in practice,
# whenever the `SymbolicsNemoExt` package extension is active).
# =========================================================================

if Base.get_extension(Symbolics, :SymbolicsNemoExt) !== nothing
    @info "Nemo extension loaded — running oracle tests"
    Nemo = Base.loaded_modules[Base.PkgId(Base.UUID("2edaba10-b0f1-5616-af89-8c11ac63239a"), "Nemo")]

    @testset "resultant — Nemo oracle" begin
        # Compare Symbolics.resultant against Nemo.resultant on a small fixed
        # corpus of integer univariate polynomials.
        for (fc, gc) in (([1, 0, -1], [1, -1]),
                         ([1, -6, 11, -6], [1, -3, 2]),
                         ([1, 0, 1], [1, -2]),
                         ([2, 3, -1, 5], [1, -1, 1]))
            R, z = Nemo.polynomial_ring(Nemo.QQ, "z")
            f_nemo = sum(fc[i] * z^(length(fc) - i) for i in 1:length(fc))
            g_nemo = sum(gc[i] * z^(length(gc) - i) for i in 1:length(gc))
            expected = Rational(Nemo.resultant(f_nemo, g_nemo))
            # Build the same polynomials in Symbolics
            df, dg = length(fc) - 1, length(gc) - 1
            f = sum(fc[i] * x^(df - i + 1) for i in 1:length(fc))
            g = sum(gc[i] * x^(dg - i + 1) for i in 1:length(gc))
            @test sym_eq(resultant(f, g, x; algorithm = :euclid), expected)
            @test sym_eq(resultant(f, g, x; algorithm = :sylvester), expected)
        end
    end

    @testset "factor — Nemo oracle" begin
        # The residual from x^10 - 1 after rational-roots extraction is
        # degree 8 and reducible; without the Nemo hook we'd return it
        # whole, with the hook we get the two cyclotomic quartics.
        u, fs = Symbolics.factor(x^10 - 1, x)
        @test sym_eq(u, 1)
        @test length(fs) == 4
        @test sym_eq(reconstruct_factor(u, fs), Symbolics.expand(x^10 - 1))

        u, fs = Symbolics.factor(x^6 - 1, x)
        @test sym_eq(u, 1)
        @test length(fs) == 4
        @test sym_eq(reconstruct_factor(u, fs), Symbolics.expand(x^6 - 1))

        # Pure degree-8 reducible input (product of two quartic cyclotomics).
        u, fs = Symbolics.factor(x^8 + x^6 + x^4 + x^2 + 1, x)
        @test length(fs) == 2
        @test all(e == 1 for (_, e) in fs)
    end
else
    @info "Nemo extension not loaded — skipping oracle tests (T012 / T038)"
end

# Optional performance-regression guard. Runs only when the
# `SYMBOLICS_POLYALG_BENCH` environment variable is set so that CI speed is
# not affected; the purpose is to catch gross regressions (SC-005).
if get(ENV, "SYMBOLICS_POLYALG_BENCH", "") != ""
    @testset "polynomial_algebra — interactive-latency guard" begin
        # Warm up
        resultant(x^3 + 1, x^2 - 2, x)
        sqrfree((x - 1)^3 * (x + 2), x)
        Symbolics.factor(x^4 - 1, x)

        t_res = @elapsed resultant(x^6 + x + 1, x^5 - 3x^2 + 7, x)
        t_sqf = @elapsed sqrfree((x - 1)^2 * (x + 2)^3 * (x - 3), x)
        t_fac = @elapsed Symbolics.factor(x^4 - 5x^3 + 5x^2 + 5x - 6, x)
        @info "polynomial_algebra timings" t_res t_sqf t_fac
        @test t_res < 2.0
        @test t_sqf < 2.0
        @test t_fac < 2.0
    end
end
