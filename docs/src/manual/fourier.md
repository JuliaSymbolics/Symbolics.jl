# Signal Functions

Symbolics provides the elementary signals used in Fourier analysis. Apart from `sinc` and
`cosc`, which come from Base, they are public but not exported, so refer to them as
`Symbolics.rect` or import them explicitly.

```@example signals
using Symbolics
using Symbolics: heaviside, rect, tri

@variables t
D = Differential(t)
expand_derivatives(D(rect(t)))
```

```@example signals
f = build_function(heaviside(t) * tri(t), t; expression = Val{false})
f.([-0.5, 0.0, 0.5])
```

| Function | Value | Derivative |
|:---------|:------|:-----------|
| `heaviside(t)` | `0`, `1/2` at `t = 0`, `1` | `dirac(t)` |
| `rect(t)` | `1` for `abs(t) < 1/2`, `1/2` at `abs(t) = 1/2`, else `0` | `dirac(t + 1/2) - dirac(t - 1/2)` |
| `tri(t)` | `max(1 - abs(t), 0)` | `rect(t + 1/2) - rect(t - 1/2)` |
| `dirac(t)`, `dirac(t, n)` | `0` for `t != 0`, `DomainError` at `0` | `dirac(t, n + 1)` |
| `diraccomb(t, T)` | `0` away from multiples of `T`, `DomainError` on them | left unevaluated |
| `sinc(t)` | `sin(π t) / (π t)`, `1` at `t = 0` | `cosc(t)` |

`sinc` is Julia's normalized `sinc`. The second derivative of `sinc` is expressed as
`-π^2 * sinc(t) - 2 * cosc(t) / t`, which evaluates to `NaN` at `t = 0`.

All of these except `sinc` and `cosc` are registered with
[`@register_discontinuity`](@ref Symbolics.@register_discontinuity).

```@docs
Symbolics.heaviside
Symbolics.dirac
Symbolics.rect
Symbolics.tri
Symbolics.diraccomb
```
