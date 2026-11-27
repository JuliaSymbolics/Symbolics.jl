# Fourier Transforms

[`fourier`](@ref) and [`inverse_fourier`](@ref) compute symbolic Fourier transforms. The
transform is an unevaluated operator, [`FourierTransform`](@ref), in the same way that
[`Differential`](@ref) is: applying it builds a term, and [`expand_transforms`](@ref) applies
the transform rules to it, as [`expand_derivatives`](@ref) applies the derivative rules.

```julia
using Symbolics
@variables t ω t₀
@syms f(t)

fourier(3f(t), t => ω)       # 3F(ω)
fourier(f(t - t₀), t => ω)   # e^{-iωt₀} F(ω)
fourier(f(-2t), t => ω)      # F(-ω/2)/2
fourier(f(t)^2, t => ω)      # stays unevaluated
```

The rules are linearity and the shift, scaling and reversal of a symbolic function. Anything
else stays an exact, unevaluated transform term. Results are `Complex{Num}` in Cartesian form,
so a transform term appears through `real(FourierTransform(t => ω)(f(t)))` and
`imag(FourierTransform(t => ω)(f(t)))`. The `convention` keyword selects `:angular` (the
default), `:ordinary` or `:unitary`.

```@docs
FourierTransform
Symbolics.InverseFourierTransform
fourier
inverse_fourier
expand_transforms
```
