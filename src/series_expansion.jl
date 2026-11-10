const SeriesCoeff = Union{Int, Rational{Int}, BigInt, Rational{BigInt}, Float64, SymbolicT}
const SeriesDict = Dict{Rational{Int}, SeriesCoeff}

"""
    TruncatedSeries(coeffs::SeriesDict, order::Rational{Int}[, var])

A truncated power/Laurent/Puiseux series `Σ coeffs[e] * var^e + O(var^order)`.

`coeffs` maps an exponent `e` (`Rational{Int}`) to its coefficient (an `Int`, a
`Rational{Int}`, a `Float64`, or a symbolic value); only exponents below `order` are kept and zero
coefficients are dropped. `var` is only used for printing: the expansion variable, or
`x - x0` for an expansion about `x0` (default `x`). Immutable: arithmetic (`+`, `*`)
returns new values.
"""
struct TruncatedSeries
    coeffs::SeriesDict
    order::Rational{Int}
    var::String
end

TruncatedSeries(coeffs::SeriesDict, order::Rational{Int}) = TruncatedSeries(coeffs, order, "x")

_ts_coeff(c::Int) = c
_ts_coeff(c::Rational{Int}) = c
_ts_coeff(c::Float64) = c
function _ts_coeff(c)
    u = value(c)
    u isa Union{Integer, Rational} && return _narrow_coeff(u)
    u isa AbstractFloat && return Float64(u)
    return u::SeriesCoeff
end

const _Exact = Union{Integer, Rational}
const _Real = Union{Integer, Rational, AbstractFloat}

_has_float(c::Number) = c isa AbstractFloat
function _has_float(c)
    u = unwrap_const(unwrap(c))
    u isa Number && return u isa AbstractFloat
    iscall(u) || return false
    return any(_has_float, arguments(u))
end

"""
    _float_like(a, b)

`a` for arithmetic with `b`. An exact `a` too large for `Int` becomes `Float64` when `b`
holds a float, so float inputs never promote to `BigFloat`.
"""
function _float_like(a, b)
    a isa _Exact || return a
    n = _narrow_coeff(a)
    n isa Union{Int, Rational{Int}} && return n
    return _has_float(b) ? Float64(n) : n
end

_ts_div(a::_Exact, b::_Exact) = a // b
_ts_div(a::AbstractFloat, b::_Real) = Float64(a) / Float64(b)
_ts_div(a::_Exact, b::AbstractFloat) = Float64(a) / Float64(b)
_ts_div(a, b) = _float_like(a, b) / _float_like(b, a)

_ts_mul(a::_Exact, b::_Exact) = (big(a) * big(b))
_ts_mul(a::AbstractFloat, b::_Real) = Float64(a) * Float64(b)
_ts_mul(a::_Exact, b::AbstractFloat) = Float64(a) * Float64(b)
_ts_mul(a, b) = _float_like(a, b) * _float_like(b, a)

_ts_add(a::_Exact, b::_Exact) = (big(a) + big(b))
_ts_add(a::AbstractFloat, b::_Real) = Float64(a) + Float64(b)
_ts_add(a::_Exact, b::AbstractFloat) = Float64(a) + Float64(b)
_ts_add(a, b) = _float_like(a, b) + _float_like(b, a)

function _ts_clean(d, order::Rational{Int})
    out = SeriesDict()
    for (e, c) in d
        e >= order && continue
        _iszero(c) && continue
        out[e] = _ts_coeff(c)
    end
    return out
end

_series_zero(order) = TruncatedSeries(SeriesDict(), Rational{Int}(order))
_series_one(order) = _series_monomial(1, 0, order)

function _series_shift(s::TruncatedSeries, de::Rational{Int}, order)
    o = Rational{Int}(order)
    d = SeriesDict(e + de => c for (e, c) in s.coeffs)
    return TruncatedSeries(_ts_clean(d, o), o)
end

_series_truncate(s::TruncatedSeries, order) = _series_shift(s, 0//1, min(s.order, Rational{Int}(order)))

function _series_monomial(c, e, order)
    o = Rational{Int}(order)
    er = Rational{Int}(e)
    d = (er < o && !_iszero(c)) ? SeriesDict(er => _ts_coeff(c)) : SeriesDict()
    return TruncatedSeries(d, o)
end

_series_valuation(s::TruncatedSeries) = isempty(s.coeffs) ? s.order : minimum(keys(s.coeffs))

Base.isempty(s::TruncatedSeries) = isempty(s.coeffs)

_series_coeff_at(s::TruncatedSeries, e) = get(s.coeffs, Rational{Int}(e), 0)

function Base.:+(a::TruncatedSeries, b::TruncatedSeries)
    o = min(a.order, b.order)
    d = SeriesDict()
    for (e, c) in a.coeffs
        e < o && (d[e] = c)
    end
    for (e, c) in b.coeffs
        e < o || continue
        d[e] = haskey(d, e) ? _ts_add(d[e], c) : c
    end
    return TruncatedSeries(_ts_clean(d, o), o)
end

function Base.:*(a::TruncatedSeries, b::TruncatedSeries)
    va = _series_valuation(a)
    vb = _series_valuation(b)
    o = min(a.order + vb, b.order + va)
    d = SeriesDict()
    for (e1, c1) in a.coeffs, (e2, c2) in b.coeffs
        e = e1 + e2
        e < o || continue
        prod = _ts_mul(c1, c2)
        d[e] = haskey(d, e) ? _ts_add(d[e], prod) : prod
    end
    return TruncatedSeries(_ts_clean(d, o), o)
end

function _series_scale(c, s::TruncatedSeries)
    _iszero(c) && return _series_zero(s.order)
    d = SeriesDict(e => _ts_mul(c, cc) for (e, cc) in s.coeffs)
    return TruncatedSeries(_ts_clean(d, s.order), s.order)
end

function _unit_tail(s::TruncatedSeries, e0::Rational{Int}, c0, rel::Rational{Int})
    g = SeriesDict()
    for (e, c) in s.coeffs
        e > e0 && (g[e - e0] = _ts_div(c, c0))
    end
    return TruncatedSeries(_ts_clean(g, rel), rel)
end

function _compose_terms(delta::TruncatedSeries, order)
    isempty(delta) && return 0
    v = _series_valuation(delta)
    v > 0 || error("composition needs a vanishing inner series, got valuation $v")
    return Int(ceil(order / v))
end

"""
    _compose(c, delta, order)

`Σ_k c(k) * delta^k` truncated at `order`, for a series `delta` that vanishes at 0.
"""
function _compose(c, delta::TruncatedSeries, order::Rational{Int})
    acc = _series_monomial(c(0), 0, order)
    dpow = _series_one(order)
    for k in 1:_compose_terms(delta, order)
        dpow = dpow * delta
        ck = c(k)
        _iszero(ck) || (acc = acc + _series_scale(ck, dpow))
    end
    return acc
end

function _split_const(f, arg, x, order::Rational{Int})
    e0 = _lead_exp(arg, x)
    e0 < 0 && throw(ArgumentError(
        "$f of an expression with a pole (leading x^$e0) has no Laurent/Puiseux series " *
        "at $x = 0"))
    s = _nseries(arg, x, order)
    delta = TruncatedSeries(filter(p -> !iszero(p.first), s.coeffs), s.order)
    return _series_coeff_at(s, 0), delta
end

_exp_coeff(k) = 1 // factorial(big(k))
_sin_coeff(k) = isodd(k) ? (iseven(k ÷ 2) ? 1 : -1) // factorial(big(k)) : 0
_cos_coeff(k) = iseven(k) ? (iseven(k ÷ 2) ? 1 : -1) // factorial(big(k)) : 0
_log1p_coeff(k) = isodd(k) ? 1 // k : -1 // k
_binom_coeff(q, k) = prod(q - j for j in 0:(k - 1); init = big(1) // 1) // factorial(big(k))

_fmt_exp(e::Rational{Int}) = denominator(e) == 1 ? string(numerator(e)) : string("(", e, ")")

_fmt_pow(v::String, e::Rational{Int}) = e == 1 ? v : string(v, "^", _fmt_exp(e))

_coeff_is(c, v) = !(value(c) isa AbstractFloat) && ((c isa Number && c == v) || isequal(value(c), v))

function _fmt_term(c, e::Rational{Int}, v::String)
    cc = (c isa Rational && denominator(c) == 1) ? numerator(c) : c
    if e == 0
        return string(cc)
    end
    xpart = _fmt_pow(v, e)
    if _coeff_is(cc, 1)
        return xpart
    elseif _coeff_is(cc, -1)
        return string("-", xpart)
    end
    cstr = _fmt_coeff(cc)
    return string(cstr, "*", xpart)
end

function _fmt_coeff(c)
    c isa Integer && return string(c)
    c isa Rational && return string("(", c, ")")
    u = value(c)
    if u isa SymbolicT && iscall(u)
        return string("(", c, ")")
    end
    return string(c)
end

function Base.show(io::IO, s::TruncatedSeries)
    oterm = string("O(", _fmt_pow(s.var, s.order), ")")
    if isempty(s.coeffs)
        print(io, oterm)
        return
    end
    terms = sort(collect(keys(s.coeffs)))
    io2 = IOBuffer()
    first = true
    for e in terms
        t = _fmt_term(s.coeffs[e], e, s.var)
        if first
            print(io2, t); first = false
        elseif startswith(t, "-")
            print(io2, " - ", t[2:end])
        else
            print(io2, " + ", t)
        end
    end
    print(io2, " + ", oterm)
    print(io, String(take!(io2)))
end

_series_has_var(expr, x) = SymbolicUtils.query(isequal(unwrap(x)), unwrap(expr))

_nseries(expr, x, order::Integer) = _nseries(expr, x, Rational{Int}(order))

function _nseries(expr, x, order::Rational{Int})
    u = unwrap(expr)
    ux = unwrap(x)

    isequal(u, ux) && return _series_monomial(1, 1, order)

    if !_series_has_var(u, ux)
        cval = value(u)
        if cval isa Complex && !isreal(cval)
            throw(ArgumentError(
                "series expansion is only supported for real expressions, but found a " *
                "complex constant $cval in the expression"))
        end
        return _series_monomial(cval, 0, order)
    end

    if SymbolicUtils.isadd(u)
        acc = _series_zero(order)
        for arg in arguments(u)
            acc = acc + _nseries(arg, x, order)
        end
        return acc

    elseif SymbolicUtils.ismul(u)
        args = arguments(u)
        lead = [_lead_exp(a, x) for a in args]
        n0 = sum(lead)
        n0 >= order && return _series_zero(order)
        acc = _series_one(order - n0)
        for (j, arg) in enumerate(args)
            fj = _nseries(arg, x, order - n0 + lead[j])
            acc = acc * fj
        end
        return _series_truncate(acc, order)

    elseif SymbolicUtils.isdiv(u)
        num, den = arguments(u)
        return _nseries_div(num, den, x, order)

    elseif SymbolicUtils.ispow(u)
        base, expo = arguments(u)
        e = value(expo)
        q = _exp_rational(e)
        q === nothing && return _nseries(exp(expo * log(base)), x, order)
        s = _nseries_pow(base, q, x, order)
        return e isa AbstractFloat ? _series_float(s) : s

    elseif SymbolicUtils.isterm(u)
        op = operation(u)
        if op === complex
            throw(ArgumentError(
                "series expansion is only supported for real expressions, but got a " *
                "complex-valued expression: $(u)"))
        end
        args = arguments(u)
        length(args) == 1 || return _nseries_diff(expr, x, order)
        a = only(args)
        if op === sqrt
            return _nseries_pow(a, 1//2, x, order)
        elseif op === exp
            return _nseries_exp(a, x, order)
        elseif op === sin || op === cos
            return _nseries_trig(op, a, x, order)
        elseif op === log
            return _nseries_log(a, x, order)
        else
            return _nseries_unary(op, a, x, order)
        end
    end

    return _nseries_diff(expr, x, order)
end

"""
    _exp_rational(e)

The exponent `e` as a `Rational{Int}`, or `nothing` if it is not a real number. Float
exponents are converted with `rationalize`, so `1.5` becomes `3//2` and `0.1` becomes
`1//10`.
"""
_exp_rational(e::_Exact) = Rational{Int}(e)
_exp_rational(e::AbstractFloat) = isfinite(e) ? rationalize(Int, e) : nothing
_exp_rational(e) = nothing

_float_coeff(c) = c isa _Exact ? Float64(c) : c
_series_float(s::TruncatedSeries) =
    TruncatedSeries(SeriesDict(e => _float_coeff(c) for (e, c) in s.coeffs), s.order)

function _lead_exp(expr, x; start = 4//1, cap = 256//1)
    u = unwrap(expr); ux = unwrap(x)
    isequal(u, ux) && return 1//1
    _series_has_var(u, ux) || return 0//1
    if SymbolicUtils.ismul(u)
        return sum(_lead_exp(a, x; start, cap) for a in arguments(u))
    elseif SymbolicUtils.isdiv(u)
        num, den = arguments(u)
        return _lead_exp(num, x; start, cap) - _lead_exp(den, x; start, cap)
    elseif SymbolicUtils.ispow(u)
        base, expo = arguments(u)
        q = _exp_rational(value(expo))
        q === nothing || return q * _lead_exp(base, x; start, cap)
    elseif SymbolicUtils.isterm(u) && operation(u) === sqrt
        return (1//2) * _lead_exp(only(arguments(u)), x; start, cap)
    end
    probe = Rational{Int}(start)
    while probe <= cap
        s = _nseries(expr, x, probe)
        isempty(s) || return _series_valuation(s)
        probe *= 2
    end
    return Rational{Int}(cap)
end

function _nseries_div(num, den, x, order::Rational{Int})
    dv = _lead_exp(den, x)
    nv = _lead_exp(num, x)
    nv - dv >= order && return _series_zero(order)
    recip = _nseries_reciprocal(den, x, order - nv, dv)
    num_s = _nseries(num, x, order + dv)
    return _series_truncate(num_s * recip, order)
end

function _nseries_reciprocal(den, x, order::Rational{Int}, dv::Rational{Int})
    rel_order = order + dv
    rel_order <= 0 && return _series_zero(order)
    den_s = _nseries(den, x, dv + rel_order)
    u0 = _series_coeff_at(den_s, dv)
    _iszero(u0) && error("reciprocal: denominator leading coefficient is zero (bad valuation $dv)")

    uexps = sort!(Rational{Int}[e - dv for e in keys(den_s.coeffs) if e - dv > 0])
    D = isempty(uexps) ? 1 : foldl(lcm, (denominator(e) for e in uexps); init = 1)

    r = SeriesDict()
    r[0//1] = _ts_div(1, u0)
    nmax = Int(ceil(rel_order * D)) - 1
    for n in 1:nmax
        e = (n // D)::Rational{Int}
        e >= rel_order && break
        s = 0
        for f in uexps
            f > e && break
            rem = e - f
            haskey(r, rem) || continue
            s = _ts_add(s, _ts_mul(_series_coeff_at(den_s, dv + f), r[rem]))
        end
        _iszero(s) || (r[e] = _ts_coeff(-_ts_div(s, u0)))
    end

    return _series_shift(TruncatedSeries(r, rel_order), -dv, order)
end

function _nseries_pow(base, q::Rational{Int}, x, order::Rational{Int})
    if q == 0
        return _series_one(order)
    end
    if q > 0 && denominator(q) == 1
        p = Int(q)
        v = _lead_exp(base, x)
        b = _nseries(base, x, max(order - (p - 1) * v, v + 1))
        res = nothing
        while p > 0
            isodd(p) && (res = res === nothing ? b : res * b)
            p >>= 1
            p > 0 && (b = b * b)
        end
        return _series_truncate(res, order)
    end
    e0 = _lead_exp(base, x)
    rel = order - e0 * q
    rel <= 0 && return _series_zero(order)
    base_s = _nseries(base, x, e0 + rel)
    c0 = _series_coeff_at(base_s, e0)
    _iszero(c0) && error("_nseries_pow: base leading coefficient zero")
    g = _unit_tail(base_s, e0, c0, rel)
    acc = _compose(k -> _binom_coeff(q, k), g, rel)
    return _series_shift(_series_scale(_pow_const(c0, q), acc), e0 * q, order)
end

"""
    _iroot(n, k)

The non-negative integer `r` with `r^k == n`, or `nothing` if `n` is not a perfect
`k`-th power.
"""
function _iroot(n::BigInt, k::Int)
    n < 0 && return nothing
    r = setprecision(BigFloat, 64 + 2 * ndigits(n; base = 2)) do
        round(BigInt, BigFloat(n)^(1 / BigFloat(k)))
    end
    for c in (r - 1, r, r + 1)
        c >= 0 && c^k == n && return c
    end
    return nothing
end

_root_term(v, q) = q == 1//2 ? Symbolics.term(sqrt, v) : Symbolics.term(^, v, q)

"""
    _pow_const(c0, q)

`c0^q` for a nonzero constant `c0`: exact when `c0` is rational and the root is
rational, a `Float64` when `c0` is a float, symbolic otherwise. Fractional powers of
negative numbers stay symbolic.
"""
function _pow_const(c0, q::Rational{Int})
    vc = value(c0)
    if vc isa AbstractFloat
        (vc > 0 || denominator(q) == 1) && return Float64(vc)^Float64(q)
        return _root_term(vc, q)
    end
    vc isa Union{Integer, Rational} || return _root_term(vc, q)
    b = Rational{BigInt}(vc)
    denominator(q) == 1 && return _narrow_coeff(b^numerator(q))
    rn = _iroot(numerator(b), denominator(q))
    rd = _iroot(denominator(b), denominator(q))
    (rn === nothing || rd === nothing) && return _root_term(vc, q)
    return _narrow_coeff((rn // rd)^numerator(q))
end

function _nseries_exp(arg, x, order::Rational{Int})
    c0, delta = _split_const(exp, arg, x, order)
    s = _compose(_exp_coeff, delta, order)
    return _iszero(c0) ? s : _series_scale(_sym_call(exp, c0), s)
end

function _sym_call(op, c)
    sc = op(wrap(c))
    return _fold_exact_numeric(unwrap(sc))
end

function _nseries_trig(op, arg, x, order::Rational{Int})
    c0, delta = _split_const(op, arg, x, order)
    if _iszero(c0)
        return _compose(op === sin ? _sin_coeff : _cos_coeff, delta, order)
    end
    sd = _compose(_sin_coeff, delta, order)
    cd = _compose(_cos_coeff, delta, order)
    s0, k0 = _sym_call(sin, c0), _sym_call(cos, c0)
    return op === sin ? _series_scale(s0, cd) + _series_scale(k0, sd) :
                        _series_scale(k0, cd) + _series_scale(-s0, sd)
end

function _nseries_log(arg, x, order::Rational{Int})
    e0 = _lead_exp(arg, x)
    e0 == 0 || throw(ArgumentError(
        "log of an expression with leading x^$e0 at $x = 0 has a log(x) term and is not a " *
        "Laurent/Puiseux series"))
    arg_s = _nseries(arg, x, order)
    c0 = _series_coeff_at(arg_s, 0)
    _iszero(c0) && error("_nseries_log: zero constant term")
    g = _unit_tail(arg_s, 0//1, c0, order)
    l0 = _iszero(c0 - 1) ? 0 : _sym_call(log, c0)
    return _compose(k -> k == 0 ? l0 : _log1p_coeff(k), g, order)
end

function _nseries_diff(expr, x, order::Rational{Int})
    maxord = Int(floor(order)) - 1
    maxord < 0 && return _series_zero(order)
    cs = _taylor_coeffs(unwrap(expr), unwrap(x), 0, maxord)
    foreach(c -> _check_finite(c, expr, x), cs)
    d = SeriesDict(Rational{Int}(k - 1) => c for (k, c) in enumerate(cs))
    return TruncatedSeries(_ts_clean(d, order), order)
end

"""
    _taylor_coeffs(f, y, c0, K)

Exact Taylor coefficients `[f(c0), f'(c0), f''(c0)/2!, …]` of `f` in `y` about `y = c0`,
through order `K`. Values are kept symbolic unless exact (see `_fold_exact_numeric`).
"""
# Not `taylor_coeff` from taylor.jl, for now: it only expands about 0, its default
# `fold`/`rationalize` turn irrational values into approximate rationals (#1299, e.g.
# atan(1) -> 101534659//129277943), with `fold = Val(false)` it leaves `cosh(0)/6`
# unfolded, and its `factorial(n)` overflows for n > 20.
function _taylor_coeffs(f, y, c0, K::Int)
    D = Differential(y)
    cur = expand_derivatives(f)
    out = SeriesCoeff[]
    for k in 0:K
        raw = value(substitute_in_deriv(cur, y => c0; fold = Val(false)))
        push!(out, _ts_coeff(_fold_exact_numeric((1 // factorial(big(k))) * raw)))
        k < K && (cur = expand_derivatives(D(cur)))
    end
    return out
end

_is_finite_coeff(c::Number) = isfinite(c)
function _is_finite_coeff(c)
    u = unwrap(c)
    u isa Number && return _is_finite_coeff(u)
    SymbolicUtils.isconst(u) && return _is_finite_coeff(unwrap_const(u))
    iscall(u) || return true
    args = arguments(u)
    all(_is_finite_coeff, args) || return false
    if all(a -> unwrap_const(unwrap(a)) isa Number, args)
        r = try
            operation(u)(map(a -> unwrap_const(unwrap(a)), args)...)
        catch
            return false
        end
        return !(r isa Number) || _is_finite_coeff(r)
    end
    return true
end

function _check_finite(c, expr, x)
    _is_finite_coeff(c) && return nothing
    throw(ArgumentError(
        "$expr is not analytic at $x = 0 (a Taylor coefficient is $c); its series " *
        "expansion is not supported"))
end

const _NONSMOOTH_AT_ZERO = (abs, sign, floor, ceil, round, trunc)

"""
    _pole_rewrites

Functions with a pole where their argument vanishes, mapped to an equivalent quotient of
functions that are analytic there. Used only when the Taylor coefficients of the function
are not finite, so results at other points are unchanged.
"""
const _pole_rewrites = Dict{Any, Function}(
    cot => y -> cos(y) / sin(y),
    csc => y -> 1 / sin(y),
    coth => y -> cosh(y) / sinh(y),
    csch => y -> 1 / sinh(y),
)

"""
    _real_roots

Real `k`-th root functions (odd `k`), mapped to `k`. Expanded as
`sign(c) * (sign(c) * a)^(1//k)` when the argument `a` vanishes or has a pole, where `c`
is its leading coefficient.
"""
const _real_roots = Dict{Any, Int}(cbrt => 3)

function _nseries_real_root(a, k::Int, x, order::Rational{Int})
    e0 = _lead_exp(a, x)
    c = value(_series_coeff_at(_nseries(a, x, e0 + 1), e0))
    c isa Real && c < 0 && return _series_scale(-1, _nseries_pow(-a, 1 // k, x, order))
    return _nseries_pow(a, 1 // k, x, order)
end

function _nseries_unary(g, arg, x, order::Rational{Int})
    k = get(_real_roots, g, nothing)
    k === nothing || _lead_exp(arg, x) == 0 || return _nseries_real_root(arg, k, x, order)
    c0, delta = _split_const(g, arg, x, order)
    g in _NONSMOOTH_AT_ZERO && _iszero(c0) && throw(ArgumentError(
        "$(g)($(arg)) is not smooth at $x = 0 and has no series expansion there"))
    y = Symbolics.variable(:__series_y)
    cs = _taylor_coeffs(g(y), unwrap(y), c0, _compose_terms(delta, order))
    if !all(_is_finite_coeff, cs)
        r = get(_pole_rewrites, g, nothing)
        r === nothing || return _nseries(r(arg), x, order)
        foreach(c -> _check_finite(c, "$(g)($(arg))", x), cs)
    end
    return _compose(k -> cs[k + 1], delta, order)
end

_fits_int(n) = typemin(Int) <= n <= typemax(Int)
function _narrow_coeff(c)
    if c isa Integer
        return _fits_int(c) ? Int(c) : c
    elseif c isa Rational
        n, d = numerator(c), denominator(c)
        d == 1 && return _fits_int(n) ? Int(n) : big(n)
        return (_fits_int(n) && _fits_int(d)) ? Rational{Int}(c) : c
    end
    return c
end

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
        elseif r isa AbstractFloat && isfinite(r) && any(a -> a isa AbstractFloat, nums)
            return Float64(r)
        elseif r isa AbstractFloat && isfinite(r) && isinteger(r)
            return Integer(r)
        end
    end
    sym_args = map(a -> a isa Number ? SConst(a) : a, folded)
    return maketerm(typeof(u), op, sym_args, metadata(u))
end

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
symbolic rather than being approximated as floats). Float inputs give `Float64`
coefficients wherever a float enters the computation; a float exponent is converted with
`rationalize` (`x^0.5` is expanded as `x^(1//2)`). Complex-valued expressions are not
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
x^(1//2) + O(x^2)

julia> expand_series(sin(x), x, 6)        # ordinary Taylor series
x + (-1//6)*x^3 + (1//120)*x^5 + O(x^6)
```
"""
function expand_series(expr, x, order; x0 = 0)
    e = unwrap(expr)
    xu = unwrap(x)
    if !iszero(x0)
        e = substitute_in_deriv(e, xu => xu + unwrap(x0))
    end
    s = _nseries(e, x, Rational{Int}(order))
    return TruncatedSeries(s.coeffs, s.order, _shifted_var(x, x0))
end

function _shifted_var(x, x0)
    iszero(x0) && return string(x)
    c = value(x0)
    c isa Real && c < 0 && return string("(", x, " + ", -c, ")")
    return string("(", x, " - ", c isa Number ? c : _fmt_coeff(c), ")")
end
