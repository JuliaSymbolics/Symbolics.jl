# Signal primitives. They are `public` but not exported: names such as `rect` and `tri`
# would clash with user code.

"""
    heaviside(t)

Heaviside step function: `0` for `t < 0`, `1/2` at `t = 0` and `1` for `t > 0`.
Its derivative is [`dirac`](@ref Symbolics.dirac).
"""
heaviside(t::Real) = (1 + sign(t)) / 2

# Integers and rationals wrap or overflow, so only compare them, never do arithmetic first.
const ExactReal = Union{Integer, Rational}

"""
    rect(t)

Rectangular pulse: `1` for `|t| < 1/2`, `1/2` at `|t| = 1/2` and `0` otherwise.
"""
rect(t::Real) = heaviside(one(t) / 2 - abs(t))
function rect(t::ExactReal)
    value = one(float(t))
    -0.5 < t < 0.5 && return value
    return t == 0.5 || t == -0.5 ? value / 2 : zero(value)
end

"""
    tri(t)

Triangular pulse `max(1 - |t|, 0)`.
"""
function tri(t::Real)
    # Subtract only inside the pulse: unsigned and `typemin` arguments would wrap.
    return isnan(t) || -1 < t < 1 ? one(t) - abs(t) : zero(t)
end

# Root and side functions of the discontinuities. The `ExactReal` methods only need the
# sign of the root, and subtract only where the result cannot wrap or overflow.
_rect_root(t) = abs(t) - 1 // 2
_rect_root(t::ExactReal) = (t < -0.5 || t > 0.5) - (-0.5 < t < 0.5)
# changes sign at each of the kinks -1, 0 and 1
_tri_root(t) = t * (t^2 - 1)
_tri_root(t::ExactReal) = (t > 1) - (t < -1) + (-1 < t < 0) - (0 < t < 1)
_tri_left(t) = ifelse(t < 0, 0, 1 - t)
_tri_left(t::ExactReal) = 0 <= t < 1 ? one(t) - t : zero(t)
_tri_right(t) = ifelse(t <= 0, 1 + t, 0)
_tri_right(t::ExactReal) = -1 < t <= 0 ? one(t) + t : zero(t)

# `(sinpi(x), cospi(x))`, without rounding an exact `x` first.
_exact_sincospi(x) = (sinpi(x), cospi(x))
function _exact_sincospi(x::Integer)
    unit = one(float(x))
    return (zero(unit), isodd(x) ? -unit : unit)
end
function _exact_sincospi(x::Rational)
    r = rem(big(x), 2)  # `sinpi` rounds a rational before reducing it
    # enough bits to keep `r` apart from the zeros of the sine and cosine
    bits = precision(BigFloat) + ndigits(denominator(r); base = 2)
    b = BigFloat(r; precision = bits)
    F = typeof(float(x))
    return (F(sinpi(b)), F(cospi(b)))
end

# Second derivative of `sinc`.
function dcosc(x::Number)
    y = float(x)
    y isa Real && isinf(y) && return zero(y)
    if !(abs(y) < 0.5)  # also takes NaN, on which the series would not terminate
        s, c = _exact_sincospi(x)
        # `-π^2 * sinc(x) - 2 * cosc(x) / x`, divided step by step so nothing overflows.
        # `-π` alone would be a `Float64`, whatever the type of `s`
        return (-(π * s) - 2 * (c - s / π / y) / y) / y
    end
    # The closed form is 0/0 at the origin, so sum the Taylor series near it.
    π² = abs2(oftype(real(y), π))
    u² = π² * y^2
    power, total, previous, k = -one(u²), -one(u²) / 3, zero(u²), 1
    while total != previous
        power *= -u² / ((2k - 1) * 2k)
        k += 1
        previous, total = total, total + power / (2k + 1)
    end
    return π² * total
end

# Derivative of `tri`. Only compares, so nothing is rounded away or wraps.
function dtri(t::Real)
    isnan(t) && return float(t)
    twice = (-1 <= t < 0) + (-1 < t < 0) - (0 < t <= 1) - (0 < t < 1)
    return oftype(float(t), twice / 2)
end

# Value of an impulse located where `r == 0`.
function _impulse(r::Real, t::Real)
    iszero(r) && throw(DomainError(t, "a Dirac impulse is a distribution and has no value here"))
    return isnan(r) ? float(r) : zero(float(r))
end

"""
    dirac(t)
    dirac(t, n)

Dirac delta distribution, or its `n`-th derivative for an integer `n >= 0`.

Numerically it is `0` for `t != 0` and throws a `DomainError` at `t = 0`, so that generated
code never silently produces `Inf`.
"""
dirac(t::Real) = _impulse(t, t)
function dirac(t::Real, n::Integer)
    n >= 0 || throw(DomainError(n, "the derivative order of dirac must be non-negative"))
    return dirac(t)
end

"""
    diraccomb(t, T)

Dirac comb with period `T > 0`, the sum of `dirac(t - k * T)` over all integers `k`.

Numerically it is `0` away from the impulses and throws a `DomainError` where `t` is a
multiple of `T`. That test is exact: a floating-point `t` within rounding error of a
multiple of `T` can evaluate to `0`.
"""
function diraccomb(t::Real, T::Real)
    T > 0 || throw(DomainError(T, "the period of diraccomb must be positive"))
    return _impulse(_comb_remainder(t, T), t)
end

# Zero exactly when `t` is a multiple of `T`.
function _comb_remainder(t::Real, T::Real)
    # `rem` of integers or rationals can overflow, and rounds them when mixed with floats
    stored = Union{ExactReal, AbstractFloat}
    if (t isa ExactReal || T isa ExactReal) && t isa stored && T isa stored &&
            isfinite(t) && isfinite(T)
        F = float(promote_type(typeof(t), typeof(T)))
        return iszero(rem(Rational{BigInt}(t), Rational{BigInt}(T))) ? zero(F) : one(F)
    end
    return rem(float(t), float(T))
end

@register_symbolic heaviside(t)
@register_symbolic rect(t)
@register_symbolic tri(t)
@register_symbolic dirac(t)
# The derivative order is a literal integer, so only `t` may be symbolic.
@register_symbolic dirac(t, n::Integer) false true begin
    signatures = [(Num, Integer), (SymbolicT, Integer)]
end
@register_symbolic diraccomb(t, T)
# `sinc` of a `Num` would otherwise be traced into `sinpi(t) / (π * t)`, which is `NaN` at 0.
@register_symbolic Base.sinc(x)
@register_symbolic Base.cosc(x)
@register_symbolic dcosc(x)
@register_symbolic dtri(t)

@register_derivative heaviside(t) 1 dirac(t)
@register_derivative rect(t) 1 dirac(t + 1 // 2) - dirac(t - 1 // 2)
@register_derivative tri(t) 1 dtri(t)
@register_derivative dtri(t) 1 dirac(t + 1) - 2 * dirac(t) + dirac(t - 1)
@register_derivative dirac(t) 1 dirac(t, 1)
@register_derivative dirac(t, n) 1 dirac(t, Base.Checked.checked_add(unwrap_const(n), 1))
# The derivative of `sinc` is `cosc`, which comes from DiffRules.
# `dcosc` has no derivative rule: higher derivatives stay unevaluated.
@register_derivative Base.cosc(x) 1 dcosc(x)

# `dirac` and `diraccomb` are deliberately not registered: they are zero on both sides of an
# impulse, so rewriting them as piecewise functions would silently replace them by zero.
@register_discontinuity heaviside(t) t 0 1
@register_discontinuity rect(t) _rect_root(t) 1 0
@register_discontinuity tri(t) _tri_root(t) _tri_left(t) _tri_right(t)

@public heaviside, dirac, rect, tri, diraccomb
