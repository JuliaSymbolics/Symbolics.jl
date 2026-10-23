# Taylor Series

For a real example of how to use the Taylor series functionality, see [this tutorial](@ref perturb_alg).

```@docs
series
taylor
taylor_coeff
```

## Laurent / Puiseux series expansion

While [`taylor`](@ref) and [`taylor_coeff`](@ref) require an *analytic* (ordinary Taylor)
series and raise an error at a singularity, [`expand_series`](@ref) returns the full
truncated expansion of an expression, automatically selecting:

  * a **Taylor** series when the expression is analytic at the point (non-negative integer
    exponents),
  * a **Laurent** series when it has a **pole** (negative integer exponents), and
  * a **Puiseux** series when it has a **branch point** (fractional exponents, e.g. `sqrt`).

Coefficients are kept **exact**, including irrational ones (`sqrt(2)`, `sin(1)`, `exp(1)`
are preserved symbolically rather than approximated as floating-point numbers).

```@docs
expand_series
leadterm
```

### Example

```julia
julia> using Symbolics

julia> @variables x;

julia> expand_series(exp(x)/x, x, 3)        # Laurent: simple pole at 0
x^-1 + 1 + (1//2)*x + (1//6)*x^2 + O(x^3)

julia> expand_series(sqrt(x), x, 2)         # Puiseux: branch point at 0
x^1//2 + O(x^2)

julia> leadterm(1 - cos(x), x)              # (coefficient, exponent)
(1//2, 2//1)
```

### Limitations

`expand_series` / `leadterm` compute a single-variable series around a *point on the real
line*. The following mathematical cases are **not** supported and are either rejected with
an informative `ArgumentError` or outside the current scope. They are documented here so
the limits are explicit rather than silently wrong:

  * **Complex-valued expressions** are rejected. There is no real series in `x` for e.g.
    `exp(im*x)` or any expression carrying a non-real constant; these raise an
    `ArgumentError` rather than returning a complex expansion.
  * **Essential singularities** are rejected. A pole inside a transcendental function —
    e.g. `exp(1/x)`, `sin(1/x)`, `cos(1/x)` — has no Laurent or Puiseux series and raises
    an `ArgumentError`.
  * **Logarithmic terms in the result are not representable.** An expansion whose result
    contains a `log(x)` term — e.g. `log(x)` itself, `log(1/x)`, or any `log` of an
    expression that vanishes or blows up at the point — is rejected, because a
    power/Laurent/Puiseux series cannot carry a `log(x)` term. (`log(1 + x)`, whose
    argument is nonzero and finite at `0`, expands normally.)
  * **Branch cuts and multivalued functions** are handled only through the principal
    branch implied by the chosen fractional powers / `sqrt`. There is no selection of a
    particular branch, no tracking of branch cuts, and no directional (one-sided)
    expansion. Fractional powers of a *negative or complex* constant are therefore not
    supported.
  * **Multivariate (joint) series** are not supported. `expand_series` expands in a single
    variable; other symbols appear only inside the coefficients, not as a second expansion
    variable. There is no simultaneous `O(xⁱyʲ)` truncation.
  * **Expansion about `x = ∞`** is not supported directly. Substitute `x → 1/t` and expand
    about `t = 0` if an asymptotic expansion at infinity is required.

When an expression is analytic and these limitations do not apply, `expand_series` agrees
with the Taylor results from [`taylor_coeff`](@ref).
