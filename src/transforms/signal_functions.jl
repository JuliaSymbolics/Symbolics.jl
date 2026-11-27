# Signal primitives. They are `public` but not exported: names such as `rect` and `tri`
# would clash with user code.

"""
    heaviside(t)

Heaviside step function: `0` for `t < 0`, `1/2` at `t = 0` and `1` for `t > 0`.
Its derivative is [`dirac`](@ref Symbolics.dirac).
"""
heaviside(t::Real) = (1 + sign(t)) / 2

"""
    rect(t)

Rectangular pulse: `1` for `|t| < 1/2`, `1/2` at `|t| = 1/2` and `0` otherwise.
"""
rect(t::Real) = heaviside(one(t) / 2 - abs(t))

"""
    tri(t)

Triangular pulse `max(1 - |t|, 0)`.
"""
tri(t::Real) = max(one(t) - abs(t), zero(t))

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

Numerically it is `0` away from the impulses and throws a `DomainError` where
`rem(t, T) == 0`. In floating point that test is exact: `t` within rounding error of a
multiple of `T` can evaluate to `0`.
"""
function diraccomb(t::Real, T::Real)
    T > 0 || throw(DomainError(T, "the period of diraccomb must be positive"))
    return _impulse(rem(t, T), t)
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

@register_derivative heaviside(t) 1 dirac(t)
@register_derivative rect(t) 1 dirac(t + 1 // 2) - dirac(t - 1 // 2)
@register_derivative tri(t) 1 rect(t + 1 // 2) - rect(t - 1 // 2)
@register_derivative dirac(t) 1 dirac(t, 1)
@register_derivative dirac(t, n) 1 dirac(t, unwrap_const(n) + 1)
# The derivative of `sinc` is `cosc`, which comes from DiffRules.
@register_derivative Base.cosc(x) 1 -(π^2) * sinc(x) - 2 * cosc(x) / x

@register_discontinuity heaviside(t) t 0 1
@register_discontinuity rect(t) abs(t) - 1 // 2 1 0
# `t * (t^2 - 1)` changes sign at each of the kinks -1, 0 and 1.
@register_discontinuity tri(t) t * (t^2 - 1) ifelse(t < 0, 0, 1 - t) ifelse(t <= 0, 1 + t, 0)
@register_discontinuity dirac(t) t 0 0
@register_discontinuity diraccomb(t, T) sinpi(t / T) 0 0
# `dirac(t, n)` has its impulse where `dirac(t)` does.
(::typeof(rootfunction(dirac)))(t, n) = t
(::typeof(left_continuous_function(dirac)))(t, n) = left_continuous_function(dirac)(t)
(::typeof(right_continuous_function(dirac)))(t, n) = right_continuous_function(dirac)(t)

@public heaviside, dirac, rect, tri, diraccomb
