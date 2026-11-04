# TruncatedSeries: a sparse Laurent/Puiseux series `Σ_{e<order} coeffs[e]*x^e + O(x^order)`
# about x = 0. Exponents and `order` are `Rational{Int}`, so one type covers Taylor
# (nonneg integer exponents), Laurent (negative integers, poles) and Puiseux (fractional
# exponents, branch points). Coefficients stay exact (see `_fold_exact_numeric` for the
# #1299 rule). This file holds the type, its arithmetic, the `_nseries` expansion, and the
# public `expand_series`.

# A series coefficient is an exact Int, an exact Rational{Int}, or a symbolic value (free
# parameter or irrational constant kept symbolic, e.g. sqrt(2)). The recurrences here only
# ever produce these three, so a closed Union captures them without falling back to `Any`.
const SeriesCoeff = Union{Int, Rational{Int}, SymbolicUtils.BasicSymbolic}
const SeriesDict = Dict{Rational{Int}, SeriesCoeff}

"""
    TruncatedSeries(coeffs::SeriesDict, order::Rational{Int})

A truncated power/Laurent/Puiseux series `Σ coeffs[e] * x^e + O(x^order)` about `x = 0`.

`coeffs` maps an exponent `e` (`Rational{Int}`) to its coefficient (an `Int`, a
`Rational{Int}`, or a symbolic value); only exponents below `order` are kept and zero
coefficients are dropped. Immutable: arithmetic (`+`, `*`) returns new values.
"""
struct TruncatedSeries
    coeffs::SeriesDict
    order::Rational{Int}
end

# Coefficient zero-test for numbers and symbolics. We avoid `simplify` (expensive, and can
# crash in `simplify_fractions`); a structural `isequal(value(c), 0)` suffices for the
# bookkeeping here — genuine cancellation is handled at the Add leadterm, not on every insert.
_ts_iszero(c::Number) = iszero(c)
_ts_iszero(c) = isequal(value(c), 0)

# Normalize a coefficient to the closed `SeriesCoeff` type before storing: unwrap `Num`,
# narrow wide Integer/Rational to `Int`/`Rational{Int}`. Anything else is a bug — the final
# `::SeriesCoeff` asserts loudly rather than silently widening.
_ts_coeff(c::Int) = c
_ts_coeff(c::Rational{Int}) = c
_ts_coeff(c::Integer) = Int(c)
_ts_coeff(c::Rational) = Rational{Int}(c)
_ts_coeff(c::SymbolicUtils.BasicSymbolic) = c
function _ts_coeff(c)
    u = value(c)
    u isa SeriesCoeff && return u
    u isa Integer && return Int(u)
    u isa Rational && return Rational{Int}(u)
    return u::SeriesCoeff
end

# Rational-preserving division: keep exact rationals for Int/Rational operands
# (Julia's `/` would turn 1/2 into 0.5). Fall back to `/` for symbolic/float.
_ts_div(a::Integer, b::Integer) = a // b
_ts_div(a::Union{Integer,Rational}, b::Union{Integer,Rational}) = a // b
_ts_div(a, b) = a / b

# Rational-preserving reciprocal of an integer/rational, else symbolic inverse.
_ts_inv(b::Union{Integer,Rational}) = one(b) // b
_ts_inv(b) = 1 / b

# Build a cleaned Dict (new dict, no mutation of inputs): drop terms at or beyond
# `order` and drop zero coefficients.
function _ts_clean(d, order::Rational{Int})
    out = SeriesDict()
    for (e, c) in d
        e >= order && continue
        _ts_iszero(c) && continue
        out[e] = _ts_coeff(c)
    end
    return out
end

"""
    _truncated_series(pairs, order)

Build a `TruncatedSeries` from `exponent => coefficient` pairs (or a `Dict`), truncating
at `order`, dropping zeros, and converting exponents/order to `Rational{Int}`.
"""
function _truncated_series(pairs, order)
    o = Rational{Int}(order)
    d = SeriesDict()
    for (e, c) in pairs
        d[Rational{Int}(e)] = _ts_coeff(c)
    end
    return TruncatedSeries(_ts_clean(d, o), o)
end

_truncated_series(d::Dict, order) = TruncatedSeries(_ts_clean(d, Rational{Int}(order)), Rational{Int}(order))

# A single monomial  c * x^e  known exactly to `order`.
function _series_monomial(c, e, order)
    o = Rational{Int}(order)
    er = Rational{Int}(e)
    d = (er < o && !_ts_iszero(c)) ? SeriesDict(er => _ts_coeff(c)) : SeriesDict()
    return TruncatedSeries(d, o)
end

# Valuation: smallest exponent with a nonzero coefficient; for an empty series, `order`.
_series_valuation(s::TruncatedSeries) = isempty(s.coeffs) ? s.order : minimum(keys(s.coeffs))

Base.isempty(s::TruncatedSeries) = isempty(s.coeffs)

# Coefficient of x^e (zero if absent).
_series_coeff_at(s::TruncatedSeries, e) = get(s.coeffs, Rational{Int}(e), 0)

# ── addition: O(x^min(order1, order2)) ──
function Base.:+(a::TruncatedSeries, b::TruncatedSeries)
    o = min(a.order, b.order)
    d = SeriesDict()
    for (e, c) in a.coeffs
        e < o && (d[e] = c)
    end
    for (e, c) in b.coeffs
        e < o || continue
        d[e] = haskey(d, e) ? d[e] + c : c
    end
    return TruncatedSeries(_ts_clean(d, o), o)
end

# ── multiplication ──
# Absolute precision of the product is limited by cross terms of each factor's known
# part with the other factor's O(x^order) tail: order = min(order_a + val_b, order_b + val_a).
function Base.:*(a::TruncatedSeries, b::TruncatedSeries)
    va = _series_valuation(a)
    vb = _series_valuation(b)
    o = min(a.order + vb, b.order + va)
    d = SeriesDict()
    for (e1, c1) in a.coeffs, (e2, c2) in b.coeffs
        e = e1 + e2
        e < o || continue
        prod = c1 * c2
        d[e] = haskey(d, e) ? d[e] + prod : prod
    end
    return TruncatedSeries(_ts_clean(d, o), o)
end

# Scalar multiply (constant * series): shifts nothing, scales coefficients.
function _series_scale(c, s::TruncatedSeries)
    _ts_iszero(c) && return TruncatedSeries(SeriesDict(), s.order)
    d = SeriesDict(e => c * cc for (e, cc) in s.coeffs)
    return TruncatedSeries(_ts_clean(d, s.order), s.order)
end

# Pretty-printing for debugging/tests.
# Integer-valued rational exponents/coefficients print as plain integers; the `x^0`
# factor is dropped, `x^1` prints as `x`, a unit coefficient is omitted, and fractional
# (Puiseux) exponents keep their rational form, e.g. `x^(1//2)`.

# format an exponent for display: `2//1 -> "2"`, `-1//1 -> "-1"`, `1//2 -> "1//2"`
_fmt_exp(e::Rational{Int}) = denominator(e) == 1 ? string(numerator(e)) : string(e)

# is a coefficient the exact value `v` (works for numbers and symbolics)?
_coeff_is(c, v) = (c isa Number && c == v) || isequal(value(c), v)

# render `coeff * x^e` as a human-readable term
function _fmt_term(c, e::Rational{Int})
    # coefficient that is an integer-valued rational prints as a plain integer
    cc = (c isa Rational && denominator(c) == 1) ? numerator(c) : c
    if e == 0
        return string(cc)                      # constant term: just the coefficient
    end
    xpart = e == 1 ? "x" : string("x^", _fmt_exp(e))
    if _coeff_is(cc, 1)
        return xpart                           # 1 * x^e  ->  x^e
    elseif _coeff_is(cc, -1)
        return string("-", xpart)              # -1 * x^e -> -x^e
    end
    cstr = _fmt_coeff(cc)
    return string(cstr, "*", xpart)
end

# format a coefficient: integers bare; rationals parenthesized `(1//2)`; a single symbol
# bare; a compound symbolic expression parenthesized.
function _fmt_coeff(c)
    c isa Integer && return string(c)
    c isa Rational && return string("(", c, ")")
    u = value(c)
    # a bare symbol or simple atom prints without parens; anything that is a call
    # (a sum/product/power) is parenthesized to keep the `*x^e` attachment clear.
    if u isa SymbolicUtils.BasicSymbolic && iscall(u)
        return string("(", c, ")")
    end
    return string(c)
end

function Base.show(io::IO, s::TruncatedSeries)
    if isempty(s.coeffs)
        print(io, "O(x^", _fmt_exp(s.order), ")")
        return
    end
    terms = sort(collect(keys(s.coeffs)))
    io2 = IOBuffer()
    first = true
    for e in terms
        t = _fmt_term(s.coeffs[e], e)
        if first
            print(io2, t); first = false
        elseif startswith(t, "-")
            print(io2, " - ", t[2:end])        # join `+ -a` as ` - a`
        else
            print(io2, " + ", t)
        end
    end
    print(io2, " + O(x^", _fmt_exp(s.order), ")")
    print(io, String(take!(io2)))
end

# _nseries(expr, x, order): expand `expr` into a TruncatedSeries about x = 0, keeping terms
# with exponent < `order`. Handles Laurent (poles via isdiv) and Puiseux (fractional powers
# via ^ with rational exponent). Dispatch is on the BasicSymbolic structure:
#   issym           -> the variable: x^1
#   isconst/!has_x  -> constant: c*x^0
#   isadd           -> sum of sub-series
#   ismul           -> product of sub-series
#   ispow           -> base^exponent (integer/rational exponent)
#   isdiv           -> num/den via reciprocal of the denominator's unit
#   isterm exp/sin/cos/log/sqrt -> composition series
#   else            -> differentiation fallback (analytic only)

# does `expr` contain the expansion variable `x`?
function _series_has_var(expr, x)
    u = unwrap(expr); ux = unwrap(x)
    isequal(u, ux) && return true
    iscall(u) || return false
    return any(a -> _series_has_var(a, ux), arguments(u))
end

_nseries(expr, x, order::Integer) = _nseries(expr, x, Rational{Int}(order))

function _nseries(expr, x, order::Rational{Int})
    u = unwrap(expr)
    ux = unwrap(x)

    # the variable itself
    isequal(u, ux) && return _series_monomial(1, 1, order)

    # a constant w.r.t. x (number, parameter, or any x-free subtree)
    if !_series_has_var(u, ux)
        cval = value(u)
        # Complex constants (e.g. the `im` in exp(im*x)) have no real Laurent/Taylor
        # series in x here; reject cleanly rather than failing deeper in the expansion.
        if cval isa Complex && !isreal(cval)
            throw(ArgumentError(
                "series expansion is only supported for real expressions, but found a " *
                "complex constant $cval in the expression"))
        end
        return _series_monomial(cval, 0, order)
    end

    if SymbolicUtils.isadd(u)
        acc = TruncatedSeries(SeriesDict(), order)
        for arg in arguments(u)
            acc = acc + _nseries(arg, x, order)
        end
        return acc

    elseif SymbolicUtils.ismul(u)
        # order needed from each factor depends on the others' valuations. Compute
        # each factor's leading exponent first (cheap), then expand each to the
        # absolute order it must reach so the product is correct to `order`.
        args = arguments(u)
        lead = [_series_leading_exponent(a, x) for a in args]
        n0 = sum(lead)
        # Start from the series `1` (= 1 + O(x^order)); each `*` narrows precision correctly.
        acc = _series_monomial(1, 0, order)
        # expand factor j to order:  order - (n0 - lead[j])   (so product reaches `order`)
        for (j, arg) in enumerate(args)
            fj = _nseries(arg, x, order - n0 + lead[j])
            acc = acc * fj
        end
        return acc

    elseif SymbolicUtils.isdiv(u)
        num, den = arguments(u)
        return _nseries_div(num, den, x, order)

    elseif SymbolicUtils.ispow(u)
        base, expo = arguments(u)
        e = value(expo)
        if e isa Integer || e isa Rational
            return _nseries_pow(base, Rational{Int}(e), x, order)
        else
            # variable/ symbolic exponent: b^e = exp(e*log(b))
            return _nseries(exp(expo * log(base)), x, order)
        end

    elseif SymbolicUtils.isterm(u)
        op = operation(u)
        # A `complex(re, im)` term (e.g. from exp(im*x) -> complex(cos(x), sin(x))) has no
        # real series in x; reject cleanly instead of failing on `only(arguments)` below.
        if op === complex
            throw(ArgumentError(
                "series expansion is only supported for real expressions, but got a " *
                "complex-valued expression: $(u)"))
        end
        a = only(arguments(u))
        if op === sqrt
            return _nseries_pow(a, 1//2, x, order)
        elseif op === exp
            return _nseries_exp(a, x, order)
        elseif op === sin
            return _nseries_trig(a, x, order, :sin)
        elseif op === cos
            return _nseries_trig(a, x, order, :cos)
        elseif op === log
            return _nseries_log(a, x, order)
        else
            return _nseries_diff_fallback(expr, x, order)
        end
    end

    return _nseries_diff_fallback(expr, x, order)
end

# Leading exponent of `expr` in x, used by the Mul branch to budget precision.
# Light-weight: recurse structurally; for isdiv subtract; for pow multiply; for
# transcendental-at-0 the leading exponent is 0 (f(0) finite & nonzero) unless the
# argument vanishes, in which case fall back to a short _nseries probe.
function _series_leading_exponent(expr, x)
    u = unwrap(expr); ux = unwrap(x)
    isequal(u, ux) && return 1//1
    _series_has_var(u, ux) || return 0//1
    if SymbolicUtils.ismul(u)
        return sum(_series_leading_exponent(a, x) for a in arguments(u))
    elseif SymbolicUtils.isdiv(u)
        num, den = arguments(u)
        return _series_leading_exponent(num, x) - _series_leading_exponent(den, x)
    elseif SymbolicUtils.ispow(u)
        base, expo = arguments(u)
        e = value(expo)
        if e isa Integer || e isa Rational
            return Rational{Int}(e) * _series_leading_exponent(base, x)
        end
    elseif SymbolicUtils.isterm(u) && operation(u) === sqrt
        return (1//2) * _series_leading_exponent(only(arguments(u)), x)
    end
    # add / transcendental / unknown: probe with a doubling order until a nonzero term
    # appears. A fixed order under-probes high-order cancellations (sin(x^3)-x^3 leads at
    # x^9) and clamps pole×polynomial precision. Terminates: the fallback fires only on
    # isadd/transcendental nodes, and _nseries recurses here only on strictly smaller
    # subtrees (inner Mul factors), never on `expr` itself.
    PROBE_START = 4//1
    PROBE_MAX   = 256//1   # ≥ this ⇒ treat as vanishing to that order (leading exp = order)
    probe_order = PROBE_START
    while probe_order <= PROBE_MAX
        s = _nseries(expr, x, probe_order)
        isempty(s) || return _series_valuation(s)
        probe_order *= 2
    end
    # Exhausted the cap without a nonzero term: the expression is zero (or vanishes
    # faster than PROBE_MAX). Report the deepest probed order as the leading exponent.
    return PROBE_MAX
end

# num / den: series(num) * (1/den). 1/den via reciprocal of the denominator's unit.
function _nseries_div(num, den, x, order::Rational{Int})
    dv = _series_leading_exponent(den, x)          # denominator leading exponent
    # series of den shifted to a unit: den = x^dv * unit, unit(0) != 0.
    recip = _nseries_reciprocal(den, x, order, dv)
    num_s = _nseries(num, x, order + dv)           # need dv extra orders
    return num_s * recip
end

# reciprocal series of `den` to absolute `order`, given leading exponent `dv`.
function _nseries_reciprocal(den, x, order::Rational{Int}, dv::Rational{Int})
    # Expand den to enough orders, then invert the unit via the standard recurrence.
    rel_order = order + dv     # need unit coeffs through this many relative orders
    den_s = _nseries(den, x, dv + rel_order)
    # Keep the full `rel_order` relative precision (NOT floor(order)): truncating earlier
    # dropped `dv` terms when the denominator had a pole, e.g. 1/(x²(1+x)) lost its x¹ term.
    rel_order <= 0 && return TruncatedSeries(SeriesDict(), order)
    u0 = _series_coeff_at(den_s, dv)
    _ts_iszero(u0) && error("reciprocal: denominator leading coefficient is zero (bad valuation $dv)")

    # Relative unit exponents may live on a fractional (Puiseux) grid, e.g. {1//2, 1, …}
    # for 1+√x, so step the recurrence on the actual lcm grid 0, 1/D, 2/D, … not integers.
    uexps = sort!(Rational{Int}[e - dv for e in keys(den_s.coeffs) if e - dv > 0])
    D = isempty(uexps) ? 1 : foldl(lcm, (denominator(e) for e in uexps); init = 1)

    # r[0] = 1/u0; r[e] = -(1/u0) Σ_{0<f≤e} u[f]·r[e−f]
    r = SeriesDict()
    r[0//1] = _ts_inv(u0)
    nmax = Int(ceil(rel_order * D)) - 1
    for n in 1:nmax
        e = (n // D)::Rational{Int}
        e >= rel_order && break
        s = 0
        for f in uexps
            f > e && break
            rem = e - f
            haskey(r, rem) || continue
            s = s + _series_coeff_at(den_s, dv + f) * r[rem]
        end
        _ts_iszero(s) || (r[e] = _ts_coeff(-_ts_div(s, u0)))
    end

    # 1/den = x^-dv * Σ r[e] x^e  ->  coeff of x^(e - dv)
    d = SeriesDict()
    for (e, c) in r
        _ts_iszero(c) || (d[e - dv] = c)
    end
    return TruncatedSeries(_ts_clean(d, order), order)
end

# base^q for rational q, via c0^q * x^(e0 q) * (1+g)^q binomial series.
function _nseries_pow(base, q::Rational{Int}, x, order::Rational{Int})
    if q == 0
        return _series_monomial(1, 0, order)
    end
    if q > 0 && denominator(q) == 1
        # positive integer power: repeated multiply (keeps exactness, avoids binomial on polynomials)
        p = Int(q)
        acc = _nseries(base, x, order)
        res = acc
        for _ in 2:p
            res = res * acc
        end
        return p == 1 ? acc : res
    end
    e0 = _series_leading_exponent(base, x)
    base_s = _nseries(base, x, e0 + (order - e0 * q))
    c0 = _series_coeff_at(base_s, e0)
    _ts_iszero(c0) && error("_nseries_pow: base leading coefficient zero")
    # g = base/(c0 x^e0) - 1
    rel = order - e0 * q
    rel <= 0 && return _series_monomial(_pow_const(c0, q), e0 * q, order)
    g = SeriesDict()
    for (e, c) in base_s.coeffs
        e > e0 && (g[e - e0] = _ts_div(c, c0))
    end
    gts = TruncatedSeries(_ts_clean(g, rel), rel)
    if isempty(gts)
        return _series_monomial(_pow_const(c0, q), e0 * q, order)
    end
    vg = _series_valuation(gts)
    maxk = Int(ceil(rel / vg))
    acc = _series_monomial(1, 0, rel)
    gpow = _series_monomial(1, 0, rel)
    binom = 1 // 1
    for k in 1:maxk
        binom = binom * (q - (k - 1)) // k
        gpow = gpow * gts
        acc = acc + _series_scale(binom, gpow)
    end
    pref = _pow_const(c0, q)
    shifted = SeriesDict()
    for (e, c) in acc.coeffs
        shifted[e + e0 * q] = pref * c
    end
    return TruncatedSeries(_ts_clean(shifted, order), order)
end

# c0^q as an exact rational when possible, else SYMBOLIC. A fractional power of an
# integer/rational base (e.g. 2^(1//2)) is irrational and must stay symbolic (`sqrt(2)`),
# not collapse to a Float64 (#1299 class). Build the symbolic power and fold only if exact.
function _pow_const(c0, q::Rational{Int})
    if (c0 isa Integer || c0 isa Rational) && denominator(q) == 1
        return Rational{Int}(c0)^Int(q)                 # exact integer power
    end
    # Fractional power of a constant (irrational, e.g. 2^(1//2)). The `^` term folds to a
    # Float64 at construction under SymReal, so we must NOT build it via `^`/`maketerm`.
    # Use `sqrt` for 1//2 (stays symbolic as a Term), and `Symbolics.term(^, …)` for other
    # fractional powers (which prints/stays symbolic). Collapse perfect powers via a cheap
    # exactness check first (e.g. 4^(1//2) ⇒ 2).
    vc = value(c0)
    if vc isa Union{Integer, Rational}
        exact = try
            r = float(vc)^float(q)
            isinteger(r) ? Integer(r) : nothing
        catch
            nothing
        end
        exact === nothing || return exact
    end
    return q == 1//2 ? Symbolics.term(sqrt, vc) : Symbolics.term(^, vc, q)
end

# exp(arg): if arg(0)=0 (vanishes), series is sum arg^k/k!. If arg(0)=c finite
# nonzero, exp(arg) = exp(c) * exp(arg-c) with (arg-c) vanishing.
function _nseries_exp(arg, x, order::Rational{Int})
    e0 = _series_leading_exponent(arg, x)
    # A pole in the argument (e0 < 0) makes exp(arg) an ESSENTIAL singularity (e.g.
    # exp(1/x)); it has no Laurent/Puiseux series. Reject cleanly rather than feeding the
    # pole to the differentiation fallback (which would substitute x=0 into 1/x → 1//0).
    e0 < 0 && throw(ArgumentError(
        "exp of an expression with a pole (leading x^$e0) is an essential singularity at " *
        "$x = 0 and has no Laurent/Puiseux series expansion"))
    arg_s = _nseries(arg, x, order)
    c0 = e0 <= 0 ? _series_coeff_at(arg_s, 0) : 0
    if !_ts_iszero(c0)
        # factor out exp(c0): keep exp(c0) SYMBOLIC (exp(1) must not collapse to 2.718);
        # _sym_call builds the symbolic call and folds it only if the result is exact.
        shifted = _nseries(arg - c0, x, order)
        inner = _exp_vanishing(shifted, order)
        return _series_scale(_sym_call(exp, c0), inner)
    end
    return _exp_vanishing(arg_s, order)
end

# Apply `op` to a constant, keeping the result SYMBOLIC unless it is exact. Avoids the
# #1299-class leak where `exp(1)` / `log(2)` would evaluate to a Float64. Builds the
# symbolic call `op(wrap(c))` then runs `_fold_exact_numeric` (so `exp(0)⇒1`, `log(1)⇒0`
# collapse, while `exp(1)`, `log(2)` stay symbolic).
function _sym_call(op, c)
    sc = op(wrap(c))            # symbolic (Num) call, e.g. exp(Num(1))
    return _fold_exact_numeric(unwrap(sc))
end

# exp of a series whose valuation is > 0
function _exp_vanishing(s::TruncatedSeries, order::Rational{Int})
    v = _series_valuation(s)
    v <= 0 && !isempty(s) && error("_exp_vanishing expects a vanishing argument, got valuation $v")
    maxk = isempty(s) ? 0 : Int(ceil(order / v))
    acc = _series_monomial(1, 0, order)
    spow = _series_monomial(1, 0, order)
    fact = 1
    for k in 1:maxk
        fact *= k
        spow = spow * s
        acc = acc + _series_scale(1 // fact, spow)
    end
    return acc
end

# sin/cos of `arg`, dispatched on the argument's leading exponent e0:
#   e0 < 0  -> pole inside the function: essential singularity, rejected.
#   e0 == 0 -> argument finite & nonzero at 0 (e.g. sin(x+1), sin(1+√x)): factor
#              arg = c0 + δ with v(δ) > 0 and use the angle-addition identities so the
#              δ part expands as a vanishing series (works even when arg is non-analytic,
#              e.g. 1+√x, whose derivative blows up at 0).
#   e0 > 0  -> argument vanishes (e.g. sin(x), sin(x^2)): expand via the Taylor series
#              of sin/cos composed with the argument series.
function _nseries_trig(arg, x, order::Rational{Int}, which::Symbol)
    e0 = _series_leading_exponent(arg, x)
    # A pole in the argument (e0 < 0) makes sin/cos(arg) an ESSENTIAL singularity (e.g.
    # sin(1/x)); reject cleanly instead of feeding the pole to the fallback (which would
    # substitute x=0 into 1/x → 1//0 and emit garbage coefficients).
    e0 < 0 && throw(ArgumentError(
        "$(which) of an expression with a pole (leading x^$e0) is an essential singularity " *
        "at $x = 0 and has no Laurent/Puiseux series expansion"))
    if e0 == 0
        # arg finite & nonzero at 0. Factor arg = c0 + δ, v(δ) > 0, and use angle addition
        #   sin(c0+δ) = sin(c0)cos(δ) + cos(c0)sin(δ)
        #   cos(c0+δ) = cos(c0)cos(δ) − sin(c0)sin(δ)
        # keeping sin(c0)/cos(c0) SYMBOLIC (sin(1) must not collapse to 0.841…). The δ
        # series is vanishing, so _trig_vanishing expands it directly — this handles
        # non-analytic args like 1+√x that the differentiation fallback cannot.
        arg_s = _nseries(arg, x, order)
        c0 = _series_coeff_at(arg_s, 0)
        delta = _nseries(arg - c0, x, order)
        sin_d = _trig_vanishing(delta, order, :sin)
        cos_d = _trig_vanishing(delta, order, :cos)
        sc0 = _sym_call(sin, c0)
        cc0 = _sym_call(cos, c0)
        return which === :sin ?
            _series_scale(sc0, cos_d) + _series_scale(cc0, sin_d) :
            _series_scale(cc0, cos_d) + _series_scale(-sc0, sin_d)
    end
    return _trig_vanishing(_nseries(arg, x, order), order, which)
end

# sin/cos of a series `s` whose valuation is > 0 (δ vanishes at 0), via the Taylor series
# Σ ±s^k/k!. Shared by the e0>0 path and the angle-addition e0==0 path.
function _trig_vanishing(s::TruncatedSeries, order::Rational{Int}, which::Symbol)
    if isempty(s)
        # δ ≡ 0: sin(0)=0, cos(0)=1
        return which === :sin ? TruncatedSeries(SeriesDict(), order) : _series_monomial(1, 0, order)
    end
    v = _series_valuation(s)
    v <= 0 && error("_trig_vanishing expects a vanishing argument, got valuation $v")
    maxk = Int(ceil(order / v))
    acc = TruncatedSeries(SeriesDict(), order)
    spow = _series_monomial(1, 0, order)  # s^0
    fact = 1
    for k in 0:maxk
        k > 0 && (fact *= k; spow = spow * s)
        coeff = if which === :sin
            isodd(k) ? ((((k - 1) ÷ 2) % 2 == 0) ? 1 // fact : -1 // fact) : 0 // 1
        else
            iseven(k) ? (((k ÷ 2) % 2 == 0) ? 1 // fact : -1 // fact) : 0 // 1
        end
        iszero(coeff) || (acc = acc + _series_scale(coeff, spow))
    end
    return acc
end

# log(arg): write arg = c0 * x^e0 * (1 + g). Then log = log(c0) + e0*log(x) + Σ (-1)^(k-1) g^k/k.
# The e0*log(x) term is NOT a power of x, so if e0 != 0 the result is not a Laurent/Puiseux
# series and we reject it (see the guard below).
function _nseries_log(arg, x, order::Rational{Int})
    e0 = _series_leading_exponent(arg, x)
    # log(arg) is a Laurent/Puiseux series only when arg has a nonzero finite value at 0
    # (e0 == 0). If arg vanishes (e0 > 0, e.g. log(x)) or has a pole (e0 < 0, e.g. log(1/x)),
    # the expansion contains a log(x) term, which a power/Laurent/Puiseux series cannot
    # represent; reject cleanly.
    e0 == 0 || throw(ArgumentError(
        "log of an expression with leading x^$e0 at $x = 0 has a log(x) term and is not a " *
        "Laurent/Puiseux series"))
    arg_s = _nseries(arg, x, order)
    c0 = _series_coeff_at(arg_s, 0)
    _ts_iszero(c0) && error("_nseries_log: zero constant term")
    # g = arg/c0 - 1
    g = SeriesDict()
    for (e, c) in arg_s.coeffs
        e > 0 && (g[e] = _ts_div(c, c0))
    end
    gts = TruncatedSeries(_ts_clean(g, order), order)
    acc = TruncatedSeries(SeriesDict(), order)
    const_term = _ts_iszero(c0 - 1) ? 0 : _sym_call(log, c0)
    _ts_iszero(const_term) || (acc = acc + _series_monomial(const_term, 0, order))
    if !isempty(gts)
        vg = _series_valuation(gts)
        maxk = Int(ceil(order / vg))
        gpow = _series_monomial(1, 0, order)
        for k in 1:maxk
            gpow = gpow * gts
            sign = isodd(k) ? (1 // k) : (-1 // k)
            acc = acc + _series_scale(sign, gpow)
        end
    end
    return acc
end

# Differentiation fallback for analytic expressions the above cases don't cover
# (e.g. tan, atan, sin/cos with a non-vanishing argument like sin(x+1)).
# Only valid when expr is analytic at 0 (no pole/branch); returns nonneg-order terms.
#
# EXACTNESS (JuliaSymbolics/Symbolics.jl#1299): coefficients are computed from raw
# derivatives substituted at 0 with `fold=Val(false)`, which keeps `exp(0)`, `cos(0)`,
# `sin(1)` etc. SYMBOLIC instead of collapsing them to `Float64` (the #1299 symptom: a
# folded `cos(0)==1.0` makes `(1//2)*1.0==0.5`, losing exactness). We then divide by
# `k!` with exact `//` and apply `_fold_exact_numeric`, which evaluates a function call
# on literal-number arguments ONLY when the result is exact (integer/rational, or an
# integer-valued float such as `exp(0)==1.0` ⇒ `1`). A genuinely irrational result
# (`sin(1)`, `√2`) is left SYMBOLIC rather than approximated. So this path is exact for
# both rational and irrational coefficients — no `rationalize`, no silent approximation.
function _nseries_diff_fallback(expr, x, order::Rational{Int})
    maxord = Int(floor(order)) - 1
    maxord < 0 && return TruncatedSeries(SeriesDict(), order)
    D = Differential(x)
    cur = expand_derivatives(unwrap(expr))
    d = SeriesDict()
    for k in 0:maxord
        raw = value(substitute_in_deriv(cur, x => 0; fold = Val(false)))
        c = _fold_exact_numeric((1 // factorial(big(k))) * raw)
        c = _narrow_coeff(c)
        _ts_iszero(c) || (d[Rational{Int}(k)] = c)
        k < maxord && (cur = expand_derivatives(D(cur)))
    end
    return TruncatedSeries(_ts_clean(d, order), order)
end

# Narrow a Rational{BigInt} back to Rational{Int} / Int when it fits (factorial(big)
# keeps us overflow-safe but we prefer small types in the output coefficients).
function _narrow_coeff(c)
    if c isa Rational && typemin(Int) <= numerator(c) <= typemax(Int) &&
                         typemin(Int) <= denominator(c) <= typemax(Int)
        return Rational{Int}(c)
    end
    return c
end

# Selective exact constant-folding (a scoped, exactness-guarded variant of the
# SymbolicUtils rule `@rule (~f)(~x::is_literal_number) => (~f)(~x)` in simplify_rules.jl,
# which unconditionally evaluates to Float64). Here we evaluate `f(nums...)` only when
# the result is EXACT, so `exp(0)⇒1`, `cos(0)⇒1`, `sin(0)⇒0` collapse while `sin(1)`,
# `cos(1)`, `√2` stay symbolic. Recurses so `(1//2)*exp(0) ⇒ 1//2`. Kept local to the
# series-expansion code rather than changed globally in SymbolicUtils.
function _fold_exact_numeric(ex)
    u = unwrap(ex)
    u isa Number && return u
    SymbolicUtils.isconst(u) && return unwrap_const(u)
    iscall(u) || return u
    op = operation(u)
    folded = map(_fold_exact_numeric, arguments(u))
    nums = map(a -> a isa Number ? a :
                    (SymbolicUtils.isconst(unwrap(a)) ? unwrap_const(unwrap(a)) : a),
               folded)
    if !isempty(nums) && all(a -> a isa Number, nums)
        r = try
            op(nums...)
        catch
            nothing
        end
        if r isa Integer || r isa Rational
            return r
        elseif r isa AbstractFloat && isfinite(r) && isinteger(r)
            return Integer(r)
        end
        # irrational / inexact result: fall through, keep the call symbolic
    end
    # Rebuild symbolically. Plain-number args must be SConst-wrapped (per the contribution
    # notes: there is no auto-convert Number->BasicSymbolic), otherwise `maketerm` with a
    # concrete op and raw numeric args re-evaluates (e.g. 2^(1//2) ⇒ 1.414 float).
    sym_args = map(a -> a isa Number ? SConst(a) : a, folded)
    return maketerm(typeof(u), op, sym_args, metadata(u))
end

# Public series-expansion API.

"""
    expand_series(expr, x, order; x0 = 0)

Expand `expr` as a truncated series in `x` about `x = x0`, up to (but not including) the
exponent `order`, returning a `TruncatedSeries` `Σ cₑ·(x - x0)^e + O((x-x0)^order)`.

Unlike [`taylor`](@ref)/[`taylor_coeff`](@ref), which require an analytic (Taylor) series
and error on singularities, `expand_series` returns the full expansion:

  * a **Taylor** series when `expr` is analytic at `x0` (non-negative integer exponents),
  * a **Laurent** series when `expr` has a pole (negative integer exponents), and
  * a **Puiseux** series when `expr` has a branch point (fractional exponents, e.g. `sqrt`).

Coefficients are kept **exact**, including irrational ones (`sqrt(2)`, `sin(1)` stay
symbolic rather than being approximated as floats). Complex-valued expressions are not
supported and raise an `ArgumentError`.

Examples
========

```jldoctest
julia> using Symbolics

julia> @variables x
1-element Vector{Num}:
 x

julia> expand_series(exp(x)/x, x, 3)      # Laurent: simple pole at 0
x^-1 + 1 + (1//2)*x + (1//6)*x^2 + O(x^3)

julia> expand_series(sqrt(x), x, 2)       # Puiseux: branch point at 0
x^1//2 + O(x^2)

julia> expand_series(sin(x), x, 6)        # ordinary Taylor series
x + (-1//6)*x^3 + (1//120)*x^5 + O(x^6)
```
"""
function expand_series(expr, x, order; x0 = 0)
    e = unwrap(expr)
    xu = unwrap(x)
    if !iszero(x0)
        # expand about x0 by substituting x -> x + x0, so the result is in powers of (x - x0)
        e = substitute_in_deriv(e, xu => xu + unwrap(x0))
    end
    return _nseries(e, x, Rational{Int}(order))
end
