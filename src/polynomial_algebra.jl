# Polynomial algebra primitives for Symbolics.jl
#
# Exposes `resultant` for univariate polynomials in a Symbolics variable,
# with rational / symbolic coefficients treated as transcendentals. See
# `docs/src/manual/polynomial_algebra.md` for user-facing documentation and
# `specs/001-poly-factor-resultant/` for the design documents.

# ------------------------------------------------------------
# Foundational helpers
# ------------------------------------------------------------

# Return the coefficient vector (highest degree first) of `expr` viewed as a
# univariate polynomial in `var`. Coefficients are `Num`; other free symbols in
# `expr` are kept as symbolic coefficients.
#
# Throws ArgumentError if `expr` is not a polynomial in `var`.
function _univariate_coeffs(expr, var)
    e = unwrap(expr)
    v = unwrap(var)
    coeffs_dict, residual = polynomial_coeffs(e, [v])
    if !SymbolicUtils._iszero(residual)
        throw(ArgumentError("expression is not a polynomial in $(var): residual term $(wrap(residual))"))
    end
    deg = 0
    for mono in keys(coeffs_dict)
        d = _mono_power(mono, v)
        deg = max(deg, d)
    end
    result = Num[Num(0) for _ in 0:deg]
    for (mono, coeff) in coeffs_dict
        d = _mono_power(mono, v)
        result[deg - d + 1] = Num(_widen_numeric(coeff))
    end
    return result
end

# Promote literal Int / Rational{Int} coefficients to their BigInt counterparts
# before they enter the resultant PRS. Intermediate powers of numeric leading
# coefficients grow fast; without this, the Euclidean path overflows on random
# degree-5 integer polynomials.
function _widen_numeric(c)
    if c isa Integer
        return BigInt(c)
    elseif c isa Rational
        return Rational{BigInt}(c)
    else
        return c
    end
end

function _mono_power(mono, v)
    m = unwrap(mono)
    if SymbolicUtils._isone(m)
        return 0
    elseif isequal(m, v)
        return 1
    elseif iscall(m) && operation(m) === (^)
        base, expn = arguments(m)
        if isequal(unwrap(base), v)
            return Int(unwrap_const(expn))
        end
    end
    throw(ArgumentError("unexpected monomial $(mono) in polynomial coefficient dict"))
end

function _strip_leading_zeros(coeffs::Vector{Num})
    result = coeffs
    while !isempty(result) && SymbolicUtils._iszero(unwrap(simplify(result[1])))
        result = result[2:end]
    end
    return result
end

# ------------------------------------------------------------
# Resultant — Sylvester matrix
# ------------------------------------------------------------

function _resultant_sylvester(fcoeffs::Vector{Num}, gcoeffs::Vector{Num})
    df = length(fcoeffs) - 1
    dg = length(gcoeffs) - 1
    if dg == 0
        return gcoeffs[1]^df
    end
    if df == 0
        return fcoeffs[1]^dg
    end
    n = df + dg
    M = Num[Num(0) for _ in 1:n, _ in 1:n]
    for i in 1:dg
        for j in 1:(df + 1)
            M[i, i + j - 1] = fcoeffs[j]
        end
    end
    for i in 1:df
        for j in 1:(dg + 1)
            M[dg + i, i + j - 1] = gcoeffs[j]
        end
    end
    return det(M)
end

# ------------------------------------------------------------
# Resultant — extended Euclidean / polynomial remainder sequence
# ------------------------------------------------------------

function _resultant_euclid(fcoeffs::Vector{Num}, gcoeffs::Vector{Num})
    df = length(fcoeffs) - 1
    dg = length(gcoeffs) - 1
    if dg == 0
        return gcoeffs[1]^df
    end
    if df == 0
        return fcoeffs[1]^dg
    end
    if df < dg
        sign = isodd(df * dg) ? -1 : 1
        return sign * _resultant_euclid(gcoeffs, fcoeffs)
    end
    r = _poly_rem(fcoeffs, gcoeffs)
    r = _strip_leading_zeros(r)
    if isempty(r)
        return Num(0)
    end
    dr = length(r) - 1
    lc_g = gcoeffs[1]
    # Identity: Res(f, g) = lc(g)^(df - dr) * (-1)^(df * dg) * Res(g, r)
    # Derivation: prod_{g(β)=0} f(β) = prod r(β) (since f ≡ r mod g), then
    # combine the two "via roots of g" formulas and the Res(r, g)↔Res(g, r)
    # swap. See specs/001-poly-factor-resultant/research.md §R2.
    sign = isodd(df * dg) ? -1 : 1
    return sign * lc_g^(df - dr) * _resultant_euclid(gcoeffs, r)
end

# True Euclidean polynomial remainder: f = q * g + r with deg(r) < deg(g).
# Coefficients may be rational functions in the non-distinguished variables;
# the final resultant identity cancels any spurious denominators.
function _poly_rem(f::Vector{Num}, g::Vector{Num})
    r = copy(f)
    lc_g = g[1]
    while true
        r = _strip_leading_zeros(r)
        length(r) < length(g) && break
        lc_r = r[1]
        factor = lc_r / lc_g
        for j in 1:length(g)
            r[j] = r[j] - factor * g[j]
        end
    end
    return r
end

# ------------------------------------------------------------
# Public API
# ------------------------------------------------------------

"""
    resultant(f, g, var; algorithm=:euclid)

Compute the resultant of two univariate polynomials `f` and `g` with respect to
the variable `var`. Returns a `Num` (reducing to a number when both inputs have
purely numeric coefficients).

Two algorithms are selectable via the `algorithm` keyword:

  - `:euclid` (default) — polynomial remainder sequence via Euclidean division
    with symbolic-rational coefficients. Faster on dense inputs.
  - `:sylvester` — determinant of the Sylvester matrix. Simpler and more
    auditable; identical mathematical result.

Both algorithms return mathematically equal results on every valid input; any
observed disagreement is a defect.

The resultant vanishes iff `f` and `g` share a common root in the algebraic
closure of their coefficient field.

# Examples
```julia
julia> using Symbolics

julia> @variables x a b c;

julia> resultant(x^2 - 1, x - 1, x)
0

julia> resultant(x^2 + 1, x - 2, x)
5

julia> resultant(a*x^2 + b*x + c, 2a*x + b, x) # equals -a * (b^2 - 4ac)
```
"""
function resultant(f, g, var; algorithm::Symbol = :euclid)
    fcoeffs = _univariate_coeffs(f, var)
    gcoeffs = _univariate_coeffs(g, var)
    fcoeffs = _strip_leading_zeros(fcoeffs)
    gcoeffs = _strip_leading_zeros(gcoeffs)
    if isempty(fcoeffs) || isempty(gcoeffs)
        throw(ArgumentError("resultant: zero polynomial"))
    end
    r = if algorithm === :euclid
        _resultant_euclid(fcoeffs, gcoeffs)
    elseif algorithm === :sylvester
        _resultant_sylvester(fcoeffs, gcoeffs)
    else
        throw(ArgumentError("resultant: unknown algorithm $(algorithm); expected :euclid or :sylvester"))
    end
    return wrap(simplify(expand(r)))
end

"""
    discriminant(f, var)

Compute the discriminant of the univariate polynomial `f` with respect to
`var`. Returns a `Num` (reducing to a number when the coefficients are
numeric).

The discriminant vanishes iff `f` has a repeated root in the algebraic closure
of its coefficient field. The value is derived from the classical identity

    disc(f) = (-1)^(n*(n-1)/2) * resultant(f, f') / lc(f)

where `n` is the degree of `f` and `lc(f)` its leading coefficient.

Conventions:

  - `discriminant(a*x + b, x) == 1` (linear case)
  - `discriminant(c, x) == 0` for a non-zero constant `c`

# Examples
```julia
julia> using Symbolics

julia> @variables x a b c;

julia> discriminant(a*x^2 + b*x + c, x)     # b^2 - 4*a*c

julia> discriminant(x^3 - 3x + 2, x)         # double root at x = 1 → 0
0
```
"""
function discriminant(f, var)
    fcoeffs = _univariate_coeffs(f, var)
    fcoeffs = _strip_leading_zeros(fcoeffs)
    if isempty(fcoeffs)
        throw(ArgumentError("discriminant: zero polynomial"))
    end
    n = length(fcoeffs) - 1
    if n == 0
        return Num(0)
    end
    if n == 1
        return Num(1)
    end
    lc = fcoeffs[1]
    fprime = derivative(f, var)
    sign = iseven(n * (n - 1) ÷ 2) ? 1 : -1
    r = resultant(f, fprime, var)
    return wrap(simplify(expand(sign * r / lc)))
end

# ------------------------------------------------------------
# Coefficient-vector helpers for sqrfree / factor (Phases 5–6)
# ------------------------------------------------------------

# Derivative of a polynomial given by its coefficient vector (highest first).
# Returns an empty vector iff the input was a constant.
function _poly_derivative(coeffs::Vector{Num})
    n = length(coeffs) - 1
    n <= 0 && return Num[]
    return Num[coeffs[i] * (n - i + 1) for i in 1:n]
end

# Rebuild a Num polynomial from its coefficient vector (highest first).
function _coeffs_to_poly(coeffs::Vector{Num}, var)
    d = length(coeffs) - 1
    s = Num(0)
    for (i, c) in enumerate(coeffs)
        s += c * var^(d - i + 1)
    end
    return s
end

# Coefficient-vector subtraction; the result is right-aligned (constant term
# aligns at the back) and padded with zeros on the high-degree side if needed.
function _poly_subtract(a::Vector{Num}, b::Vector{Num})
    la, lb = length(a), length(b)
    n = max(la, lb)
    r = Num[Num(0) for _ in 1:n]
    for j in 1:la
        r[n - la + j] += a[j]
    end
    for j in 1:lb
        r[n - lb + j] -= b[j]
    end
    return r
end

# Divide with remainder: a = q * b + r, deg(r) < deg(b). Coefficients are
# `Num`; the division uses the coefficient field (symbolic rationals).
function _poly_divrem(a::Vector{Num}, b::Vector{Num})
    b_clean = _strip_leading_zeros(b)
    isempty(b_clean) && throw(ArgumentError("polynomial division by zero"))
    a_clean = _strip_leading_zeros(copy(a))
    da, db = length(a_clean) - 1, length(b_clean) - 1
    if da < db
        return (Num[Num(0)], isempty(a_clean) ? Num[Num(0)] : a_clean)
    end
    q = Num[Num(0) for _ in 1:(da - db + 1)]
    r = copy(a_clean)
    lc_b = b_clean[1]
    while true
        r = _strip_leading_zeros(r)
        (isempty(r) || (length(r) - 1) < db) && break
        lc_r = r[1]
        factor = lc_r / lc_b
        shift = (length(r) - 1) - db
        q[length(q) - shift] = factor
        for j in 1:length(b_clean)
            r[j] = r[j] - factor * b_clean[j]
        end
    end
    return (q, isempty(r) ? Num[Num(0)] : r)
end

# Exact division: fails (throws) if the remainder is non-zero.
function _poly_div_exact(a::Vector{Num}, b::Vector{Num})
    q, r = _poly_divrem(a, b)
    r_clean = _strip_leading_zeros(r)
    if !isempty(r_clean)
        throw(ArgumentError("polynomial division is not exact; non-zero remainder"))
    end
    return q
end

# Polynomial GCD via Euclidean recursion over the coefficient field. The
# result is returned in monic normal form (leading coefficient 1) so that
# downstream operations can safely compare lengths and factor.
function _poly_gcd(a::Vector{Num}, b::Vector{Num})
    a_local = _strip_leading_zeros(copy(a))
    b_local = _strip_leading_zeros(copy(b))
    while !isempty(b_local)
        _, r = _poly_divrem(a_local, b_local)
        r = _strip_leading_zeros(r)
        a_local = b_local
        b_local = r
    end
    isempty(a_local) && return Num[Num(1)]
    lc = a_local[1]
    return Num[ai / lc for ai in a_local]
end

# ------------------------------------------------------------
# Square-free decomposition (Yun's algorithm, characteristic 0)
# ------------------------------------------------------------

"""
    sqrfree(f, var)

Compute the square-free decomposition of the univariate polynomial `f` with
respect to `var`. Returns a tuple `(unit, factors)` where `unit::Num` is the
leading coefficient / content and `factors::Vector{Tuple{Num, Int}}` is a
list of `(p_i, e_i)` pairs such that

    f == unit * prod(p_i^e_i)

the `p_i` are pairwise coprime, each `p_i` is square-free, and every
multiplicity `e_i` is ≥ 1.

Uses Yun's algorithm over the rationals. Works for polynomials whose
coefficients are rational numbers or symbolic expressions in other
variables (treated as transcendentals).

# Examples
```julia
julia> using Symbolics

julia> @variables x;

julia> u, fs = sqrfree((x - 1)^2 * (x - 2)^3, x);
julia> u
1
julia> sort(fs; by = last)
2-element Vector{Tuple{Num, Int64}}:
 (x - 1, 2)
 (x - 2, 3)
```
"""
function sqrfree(f, var)
    fcoeffs = _univariate_coeffs(f, var)
    fcoeffs = _strip_leading_zeros(fcoeffs)
    if isempty(fcoeffs)
        throw(ArgumentError("sqrfree: zero polynomial"))
    end
    n = length(fcoeffs) - 1
    lc = fcoeffs[1]
    if n == 0
        return (wrap(simplify(lc)), Tuple{Num, Int}[])
    end
    # Work in the monic normalisation so that all intermediate polynomials
    # are monic; the leading coefficient is pulled out as the unit.
    fmonic = Num[c / lc for c in fcoeffs]

    fprime = _poly_derivative(fmonic)
    c = _poly_gcd(fmonic, fprime)
    factors = Tuple{Num, Int}[]
    if length(c) == 1
        # gcd(f, f') is a non-zero constant → f is already square-free.
        p = _coeffs_to_poly(fmonic, var)
        push!(factors, (wrap(simplify(expand(p))), 1))
        return (wrap(simplify(lc)), factors)
    end

    w = _poly_div_exact(fmonic, c)
    y = _poly_div_exact(fprime, c)
    wp = _poly_derivative(w)
    z = _poly_subtract(y, wp)

    i = 1
    safety_guard = n + 5
    while length(w) > 1 && i <= safety_guard
        p_i = _poly_gcd(w, z)
        if length(p_i) > 1
            p_sym = _coeffs_to_poly(p_i, var)
            push!(factors, (wrap(simplify(expand(p_sym))), i))
        end
        # Early exit: if p_i is trivially 1 AND z is zero, further iterations
        # will not contribute new factors.
        if length(p_i) == 1
            if all(SymbolicUtils._iszero(unwrap(simplify(zi))) for zi in z)
                break
            end
        end
        w = _poly_div_exact(w, p_i)
        y = _poly_div_exact(z, p_i)
        wp = _poly_derivative(w)
        z = _poly_subtract(y, wp)
        i += 1
    end

    return (wrap(simplify(lc)), factors)
end

# ------------------------------------------------------------
# Factor over ℚ (Phase 6 — User Story 4)
# ------------------------------------------------------------

# Attempt to coerce a Num coefficient to a Rational{BigInt}. Returns `nothing`
# if the coefficient is symbolic or non-rational.
function _try_as_rational(c)
    try
        uc = unwrap(c)
        if uc isa Integer
            return Rational{BigInt}(uc)
        elseif uc isa Rational
            return Rational{BigInt}(numerator(uc), denominator(uc))
        elseif uc isa AbstractFloat || uc isa Complex{<:AbstractFloat}
            return nothing
        elseif SU.isconst(uc)
            v = unwrap_const(uc)
            if v isa Integer
                return Rational{BigInt}(v)
            elseif v isa Rational
                return Rational{BigInt}(numerator(v), denominator(v))
            end
        end
    catch
        return nothing
    end
    return nothing
end

# Detect the first coefficient that is a floating-point literal (including
# floats hidden inside a Const). Returns the offending coefficient or nothing.
function _find_float_coefficient(coeffs::Vector{Num})
    for c in coeffs
        uc = unwrap(c)
        if uc isa AbstractFloat || uc isa Complex{<:AbstractFloat}
            return c
        end
        if SU.isconst(uc)
            v = unwrap_const(uc)
            if v isa AbstractFloat || v isa Complex{<:AbstractFloat}
                return c
            end
        end
    end
    return nothing
end

# Try to extract every coefficient as Rational{BigInt}. Returns `nothing` if
# any coefficient is symbolic; that signals the caller to fall back to
# (e.g.) marking the factor unverified rather than running a rational-roots
# sweep that wouldn't make sense.
function _try_extract_rational_vector(coeffs::Vector{Num})
    out = Rational{BigInt}[]
    for c in coeffs
        r = _try_as_rational(c)
        r === nothing && return nothing
        push!(out, r)
    end
    return out
end

# Positive divisors of an integer. Naive O(n); the integers we see here come
# from rational-root enumeration, which is interactive-grade.
function _divisors(n::Integer)
    n = abs(BigInt(n))
    n == 0 && return BigInt[]
    ds = BigInt[]
    d = BigInt(1)
    while d <= n
        n % d == 0 && push!(ds, d)
        d += 1
    end
    return ds
end

# Rational-root theorem: return all rational roots of a polynomial given by
# its Rational{BigInt} coefficient vector (highest degree first).
function _rational_roots(rat_coeffs::Vector{Rational{BigInt}})
    isempty(rat_coeffs) && return Rational{BigInt}[]
    den_lcm = reduce(lcm, denominator.(rat_coeffs); init = BigInt(1))
    int_coeffs = BigInt[numerator(c) * (den_lcm ÷ denominator(c)) for c in rat_coeffs]
    lead = int_coeffs[1]
    const_term = int_coeffs[end]
    roots = Rational{BigInt}[]
    if const_term == 0
        push!(roots, Rational{BigInt}(0))
    end
    if const_term != 0
        p_divs = _divisors(const_term)
        q_divs = _divisors(lead)
        for p in p_divs, q in q_divs
            for s in (1, -1)
                r = Rational{BigInt}(s * p, q)
                _evaluate_rational_poly(rat_coeffs, r) == 0 && push!(roots, r)
            end
        end
    end
    return unique(roots)
end

function _evaluate_rational_poly(coeffs::Vector{Rational{BigInt}}, x::Rational{BigInt})
    s = Rational{BigInt}(0)
    for c in coeffs
        s = s * x + c
    end
    return s
end

# Display-quality narrowing of a rational root: collapse `n//1` to the bare
# integer and demote `BigInt` to `Int` when it fits, so that a factor built
# from an integer root prints as `x - 1` rather than `x - 1//1`. Purely a
# presentation concern — the math is identical either way.
function _narrow_root(r::Rational{BigInt})
    if denominator(r) == 1
        n = numerator(r)
        return typemin(Int) <= n <= typemax(Int) ? Int(n) : n
    end
    if typemin(Int) <= numerator(r) <= typemax(Int) &&
       typemin(Int) <= denominator(r) <= typemax(Int)
        return Rational{Int}(Int(numerator(r)), Int(denominator(r)))
    end
    return r
end

# Synthetic division: divide a rational-coefficient polynomial by (x - r).
# Assumes r is a root; the result is one degree lower and the remainder is
# zero.
function _divide_by_linear(coeffs::Vector{Rational{BigInt}}, r::Rational{BigInt})
    n = length(coeffs) - 1
    q = Rational{BigInt}[Rational{BigInt}(0) for _ in 1:n]
    q[1] = coeffs[1]
    for i in 2:n
        q[i] = coeffs[i] + q[i - 1] * r
    end
    return q
end

# Pull all rational linear factors out of a single (monic, square-free) `Num`
# polynomial. Returns `(linear_factors, residual)` where `linear_factors` is
# a list of `Num` polynomials (each of the form `var - r`) and `residual` is
# either `nothing` (polynomial fully split) or a `Num` polynomial of degree
# ≥ 1 with no rational roots.
function _pull_rational_linear_factors(q_num::Num, var)
    coeffs_num = _univariate_coeffs(q_num, var)
    coeffs_num = _strip_leading_zeros(coeffs_num)
    length(coeffs_num) <= 1 && return (Num[], nothing)

    rat = _try_extract_rational_vector(coeffs_num)
    rat === nothing && return (Num[], q_num)

    linear_factors = Num[]
    remaining = rat
    while length(remaining) > 1
        roots = _rational_roots(remaining)
        isempty(roots) && break
        progressed = false
        for r in roots
            while length(remaining) > 1 && _evaluate_rational_poly(remaining, r) == 0
                remaining = _divide_by_linear(remaining, r)
                push!(linear_factors, wrap(simplify(expand(var - _narrow_root(r)))))
                progressed = true
            end
        end
        progressed || break
    end

    if length(remaining) <= 1
        return (linear_factors, nothing)
    end
    rem_poly = _coeffs_to_poly(Num[Num(c) for c in remaining], var)
    return (linear_factors, wrap(simplify(expand(rem_poly))))
end

# Extension hook: the Nemo ext may override this to deliver a full
# factorisation of polynomials the pure-Julia core cannot split. Returns
# `nothing` by default (Nemo not loaded) or a Vector of `(Num, Int)` factor
# pairs when Nemo has factored the residual.
_factor_with_nemo(residual, var) = nothing

"""
    factor(f, var)

Compute the factorisation of the univariate polynomial `f` over the rationals.
Returns a tuple `(unit, factors)` where `unit::Num` absorbs the leading
coefficient / integer content and `factors::Vector{Tuple{Num, Int}}` is a
list of `(p_i, e_i)` pairs with each `p_i` either irreducible over ℚ (when
the pure-Julia core can prove it, i.e., linear factors and degree-2/3
residuals with no rational root) or the best available decomposition
otherwise.

Coefficients must be rational (`Integer` / `Rational`) or symbolic in other
variables. Floating-point coefficients are not supported and raise
`ArgumentError`.

# Examples
```julia
julia> using Symbolics

julia> @variables x;

julia> u, fs = Symbolics.factor(x^4 - 1, x);  # (x-1)(x+1)(x^2+1)

julia> u, fs = Symbolics.factor(2x^2 - 2, x); # 2 * (x-1)(x+1)
```

If `Nemo` is loaded, residual higher-degree factors the core cannot prove
irreducible are passed to `Nemo.factor` via the `SymbolicsNemoExt`
extension. Without `Nemo`, degree-≥ 4 residuals without rational roots are
returned as a single non-trivially-factored entry.
"""
function factor(f, var)
    fcoeffs = _univariate_coeffs(f, var)
    fcoeffs = _strip_leading_zeros(fcoeffs)
    if isempty(fcoeffs)
        throw(ArgumentError("factor: zero polynomial"))
    end
    bad = _find_float_coefficient(fcoeffs)
    if bad !== nothing
        throw(ArgumentError("factor: coefficient $(bad) not in supported domain; floating-point coefficients are not supported"))
    end

    unit, sqfree_factors = sqrfree(f, var)
    all_factors = Tuple{Num, Int}[]
    for (q, mult) in sqfree_factors
        linear, residual = _pull_rational_linear_factors(q, var)
        for lf in linear
            push!(all_factors, (lf, mult))
        end
        residual === nothing && continue

        res_coeffs = _univariate_coeffs(residual, var)
        res_coeffs = _strip_leading_zeros(res_coeffs)
        d_res = length(res_coeffs) - 1
        if d_res <= 3
            # No rational root + degree ≤ 3 ⇒ irreducible over ℚ
            push!(all_factors, (residual, mult))
            continue
        end
        # Higher degree: try the Nemo hook
        nemo_result = _factor_with_nemo(residual, var)
        if nemo_result === nothing
            push!(all_factors, (residual, mult))
        else
            for (p, e) in nemo_result
                push!(all_factors, (p, e * mult))
            end
        end
    end
    return (unit, all_factors)
end
