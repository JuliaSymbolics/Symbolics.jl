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
| `tri(t)` | `max(1 - abs(t), 0)` | `1` on `(-1, 0)`, `-1` on `(0, 1)`, `0` outside |
| `dirac(t)`, `dirac(t, n)` | `0` for `t != 0`, `DomainError` at `0` | `dirac(t, n + 1)` |
| `diraccomb(t, T)` | `0` away from multiples of `T`, `DomainError` on them | left unevaluated |
| `sinc(t)` | `sin(π t) / (π t)`, `1` at `t = 0` | `cosc(t)` |

At its kinks the derivative of `tri` takes the mean of the two sides.

`sinc` is Julia's normalized `sinc`. Its second derivative evaluates to `-π^2/3` at `t = 0`;
the third and higher derivatives are left unevaluated.

`heaviside`, `rect` and `tri` are registered with
[`@register_discontinuity`](@ref Symbolics.@register_discontinuity). `dirac` and `diraccomb`
deliberately are not: they are zero on both sides of an impulse, so a tool that rewrites a
registered function as a piecewise function would silently replace the impulse by zero.

A root function with two sides cannot express the value at a single point. A piecewise form
reconstructed from the registration therefore takes a one-sided value exactly at a jump:
`1` instead of `1/2` for `heaviside(0)`, and `0` instead of `1/2` for `rect(±1/2)`.

```@docs
Symbolics.heaviside
Symbolics.dirac
Symbolics.rect
Symbolics.tri
Symbolics.diraccomb
```
