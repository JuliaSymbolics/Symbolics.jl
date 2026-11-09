const SeriesCoeff = Union{Int, Rational{Int}, BigInt, Rational{BigInt}, SymbolicT}
const SeriesDict = Dict{Rational{Int}, SeriesCoeff}

"""
    TruncatedSeries(coeffs::SeriesDict, order::Rational{Int}[, var])

A truncated power/Laurent/Puiseux series `Σ coeffs[e] * var^e + O(var^order)`.

`coeffs` maps an exponent `e` (`Rational{Int}`) to its coefficient (an `Int`, a
`Rational{Int}`, or a symbolic value); only exponents below `order` are kept and zero
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
_ts_coeff(c::Integer) = _narrow_coeff(c)
_ts_coeff(c::Rational) = _narrow_coeff(c)
_ts_coeff(c::SymbolicT) = c
function _ts_coeff(c)
    u = value(c)
    u isa SeriesCoeff && return _narrow_coeff(u)
    u isa Integer && return _narrow_coeff(u)
    u isa Rational && return _narrow_coeff(u)
    return u::SeriesCoeff
end

_ts_div(a::Union{Integer,Rational}, b::Union{Integer,Rational}) = a // b
_ts_div(a, b) = a / b

_ts_mul(a::Union{Integer,Rational}, b::Union{Integer,Rational}) = (big(a) * big(b))
_ts_mul(a, b) = a * b

_ts_add(a::Union{Integer,Rational}, b::Union{Integer,Rational}) = (big(a) + big(b))
_ts_add(a, b) = a + b

_ts_inv(b::Union{Integer,Rational}) = one(b) // b
_ts_inv(b) = 1 / b

function _ts_clean(d, order::Rational{Int})
    out = SeriesDict()
    for (e, c) in d
        e >= order && continue
        _iszero(c) && continue
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

_series_zero(order) = TruncatedSeries(SeriesDict(), Rational{Int}(order))
_series_one(order) = _series_monomial(1, 0, order)

function _series_truncate(s::TruncatedSeries, order)
    o = min(s.order, Rational{Int}(order))
    return TruncatedSeries(_ts_clean(s.coeffs, o), o)
end

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

_fmt_exp(e::Rational{Int}) = denominator(e) == 1 ? string(numerator(e)) : string("(", e, ")")

_fmt_pow(v::String, e::Rational{Int}) = e == 1 ? v : string(v, "^", _fmt_exp(e))

_coeff_is(c, v) = (c isa Number && c == v) || isequal(value(c), v)

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
        acc = _series_one(order)
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
            return _nseries(exp(expo * log(base)), x, order)
        end

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
        elseif op === sin
            return _nseries_trig(a, x, order, :sin)
        elseif op === cos
            return _nseries_trig(a, x, order, :cos)
        elseif op === log
            return _nseries_log(a, x, order)
        else
            return _nseries_unary(op, a, x, order)
        end
    end

    return _nseries_diff(expr, x, order)
end

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
        e = value(expo)
        if e isa Integer || e isa Rational
            return Rational{Int}(e) * _lead_exp(base, x; start, cap)
        end
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
    recip = _nseries_reciprocal(den, x, order, dv)
    num_s = _nseries(num, x, order + dv)
    return num_s * recip
end

function _nseries_reciprocal(den, x, order::Rational{Int}, dv::Rational{Int})
    rel_order = order + dv
    den_s = _nseries(den, x, dv + rel_order)
    rel_order <= 0 && return _series_zero(order)
    u0 = _series_coeff_at(den_s, dv)
    _iszero(u0) && error("reciprocal: denominator leading coefficient is zero (bad valuation $dv)")

    uexps = sort!(Rational{Int}[e - dv for e in keys(den_s.coeffs) if e - dv > 0])
    D = isempty(uexps) ? 1 : foldl(lcm, (denominator(e) for e in uexps); init = 1)

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
        _iszero(s) || (r[e] = _ts_coeff(-_ts_div(s, u0)))
    end

    d = SeriesDict()
    for (e, c) in r
        _iszero(c) || (d[e - dv] = c)
    end
    return TruncatedSeries(_ts_clean(d, order), order)
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
    base_s = _nseries(base, x, e0 + (order - e0 * q))
    c0 = _series_coeff_at(base_s, e0)
    _iszero(c0) && error("_nseries_pow: base leading coefficient zero")
    rel = order - e0 * q
    rel <= 0 && return _series_monomial(_pow_const(c0, q), e0 * q, order)
    gts = _unit_tail(base_s, e0, c0, rel)
    if isempty(gts)
        return _series_monomial(_pow_const(c0, q), e0 * q, order)
    end
    vg = _series_valuation(gts)
    maxk = Int(ceil(rel / vg))
    acc = _series_one(rel)
    gpow = _series_one(rel)
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

function _pow_const(c0, q::Rational{Int})
    if (c0 isa Integer || c0 isa Rational) && denominator(q) == 1
        return Rational{Int}(c0)^Int(q)
    end
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

function _nseries_exp(arg, x, order::Rational{Int})
    e0 = _lead_exp(arg, x)
    e0 < 0 && throw(ArgumentError(
        "exp of an expression with a pole (leading x^$e0) is an essential singularity at " *
        "$x = 0 and has no Laurent/Puiseux series expansion"))
    arg_s = _nseries(arg, x, order)
    c0 = e0 <= 0 ? _series_coeff_at(arg_s, 0) : 0
    if !_iszero(c0)
        shifted = _nseries(arg - c0, x, order)
        inner = _exp_vanishing(shifted, order)
        return _series_scale(_sym_call(exp, c0), inner)
    end
    return _exp_vanishing(arg_s, order)
end

function _sym_call(op, c)
    sc = op(wrap(c))
    return _fold_exact_numeric(unwrap(sc))
end

function _exp_vanishing(s::TruncatedSeries, order::Rational{Int})
    v = _series_valuation(s)
    v <= 0 && !isempty(s) && error("_exp_vanishing expects a vanishing argument, got valuation $v")
    maxk = isempty(s) ? 0 : Int(ceil(order / v))
    acc = _series_one(order)
    spow = _series_one(order)
    fact = big(1)
    for k in 1:maxk
        fact *= k
        spow = spow * s
        acc = acc + _series_scale(_narrow_coeff(1 // fact), spow)
    end
    return acc
end

function _nseries_trig(arg, x, order::Rational{Int}, which::Symbol)
    e0 = _lead_exp(arg, x)
    e0 < 0 && throw(ArgumentError(
        "$(which) of an expression with a pole (leading x^$e0) is an essential singularity " *
        "at $x = 0 and has no Laurent/Puiseux series expansion"))
    if e0 == 0
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

function _trig_vanishing(s::TruncatedSeries, order::Rational{Int}, which::Symbol)
    if isempty(s)
        return which === :sin ? _series_zero(order) : _series_one(order)
    end
    v = _series_valuation(s)
    v <= 0 && error("_trig_vanishing expects a vanishing argument, got valuation $v")
    maxk = Int(ceil(order / v))
    acc = _series_zero(order)
    spow = _series_one(order)
    fact = big(1)
    for k in 0:maxk
        k > 0 && (fact *= k; spow = spow * s)
        coeff = if which === :sin
            isodd(k) ? ((((k - 1) ÷ 2) % 2 == 0) ? 1 // fact : -1 // fact) : big(0) // 1
        else
            iseven(k) ? (((k ÷ 2) % 2 == 0) ? 1 // fact : -1 // fact) : big(0) // 1
        end
        iszero(coeff) || (acc = acc + _series_scale(_narrow_coeff(coeff), spow))
    end
    return acc
end

function _nseries_log(arg, x, order::Rational{Int})
    e0 = _lead_exp(arg, x)
    e0 == 0 || throw(ArgumentError(
        "log of an expression with leading x^$e0 at $x = 0 has a log(x) term and is not a " *
        "Laurent/Puiseux series"))
    arg_s = _nseries(arg, x, order)
    c0 = _series_coeff_at(arg_s, 0)
    _iszero(c0) && error("_nseries_log: zero constant term")
    gts = _unit_tail(arg_s, 0//1, c0, order)
    acc = _series_zero(order)
    const_term = _iszero(c0 - 1) ? 0 : _sym_call(log, c0)
    _iszero(const_term) || (acc = acc + _series_monomial(const_term, 0, order))
    if !isempty(gts)
        vg = _series_valuation(gts)
        maxk = Int(ceil(order / vg))
        gpow = _series_one(order)
        for k in 1:maxk
            gpow = gpow * gts
            sign = isodd(k) ? (1 // k) : (-1 // k)
            acc = acc + _series_scale(sign, gpow)
        end
    end
    return acc
end

function _nseries_diff(expr, x, order::Rational{Int})
    maxord = Int(floor(order)) - 1
    maxord < 0 && return _series_zero(order)
    cs = _taylor_coeffs(unwrap(expr), unwrap(x), 0, maxord)
    d = SeriesDict()
    for (k, c) in enumerate(cs)
        _check_finite(c, expr, x)
        _iszero(c) || (d[Rational{Int}(k - 1)] = c)
    end
    return TruncatedSeries(_ts_clean(d, order), order)
end

"""
    _taylor_coeffs(f, y, c0, K)

Exact Taylor coefficients `[f(c0), f'(c0), f''(c0)/2!, …]` of `f` in `y` about `y = c0`,
through order `K`. Values are kept symbolic unless exact (see `_fold_exact_numeric`).
"""
function _taylor_coeffs(f, y, c0, K::Int)
    D = Differential(y)
    cur = expand_derivatives(f)
    out = SeriesCoeff[]
    for k in 0:K
        raw = value(substitute_in_deriv(cur, y => c0; fold = Val(false)))
        push!(out, _ts_coeff(_narrow_coeff(_fold_exact_numeric((1 // factorial(big(k))) * raw))))
        k < K && (cur = expand_derivatives(D(cur)))
    end
    return out
end

_is_finite_coeff(c::Rational) = !iszero(denominator(c))
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

function _nseries_unary(g, arg, x, order::Rational{Int})
    e0 = _lead_exp(arg, x)
    e0 < 0 && throw(ArgumentError(
        "$(g) of an expression with a pole (leading x^$e0) is not supported: it has no " *
        "Laurent/Puiseux series expansion at $x = 0"))
    arg_s = _nseries(arg, x, order)
    c0 = _series_coeff_at(arg_s, 0)
    any(f -> g === f, _NONSMOOTH_AT_ZERO) && _iszero(c0) && throw(ArgumentError(
        "$(g)($(arg)) is not smooth at $x = 0 and has no series expansion there"))
    delta = _nseries(arg - c0, x, order)
    vd = _series_valuation(delta)
    K = isempty(delta) ? 0 : Int(ceil(order / vd))
    y = Symbolics.variable(:__series_y)
    cs = _taylor_coeffs(g(y), unwrap(y), c0, K)
    foreach(c -> _check_finite(c, "$(g)($(arg))", x), cs)
    acc = _series_monomial(cs[1], 0, order)
    dpow = _series_one(order)
    for k in 1:K
        dpow = dpow * delta
        _iszero(cs[k + 1]) || (acc = acc + _series_scale(cs[k + 1], dpow))
    end
    return acc
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
