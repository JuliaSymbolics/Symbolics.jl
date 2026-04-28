# Polynomial Algebra

Symbolics.jl provides classical univariate polynomial-algebra primitives for
polynomials whose coefficients are either rational (`Integer` / `Rational`) or
symbolic in other variables (treated as transcendentals).

## Resultant

[`Symbolics.resultant`](@ref) computes the resultant of two univariate
polynomials with respect to a distinguished variable. It vanishes exactly when
the two polynomials share a common root in the algebraic closure of the
coefficient field.

Two algorithms are selectable via the `algorithm` keyword:

- `:euclid` (default) — polynomial remainder sequence via Euclidean division.
- `:sylvester` — determinant of the Sylvester matrix.

Both return mathematically equal results; any disagreement is a defect.

```@example poly
using Symbolics
@variables x a b c d
resultant(x^2 - 1, x - 1, x)           # common root x = 1 → 0
```

```@example poly
resultant(x^2 + 1, x - 2, x)           # coprime → 5
```

```@example poly
resultant(x^2 + 1, x - 2, x; algorithm = :sylvester)
```

Symbolic coefficients are supported:

```@example poly
resultant(a*x + b, c*x + d, x)         # determinant: a*d - b*c
```

```@docs
Symbolics.resultant
```

## Discriminant

[`Symbolics.discriminant`](@ref) computes the discriminant of a univariate
polynomial. It vanishes exactly when the polynomial has a repeated root in
the algebraic closure of the coefficient field.

It is derived from the resultant via the classical identity

```math
\mathrm{disc}(f) = \frac{(-1)^{n(n-1)/2}}{\mathrm{lc}(f)} \,\mathrm{Res}(f, f'),
```

where `n` is the degree of `f`.

```@example poly
discriminant(a*x^2 + b*x + c, x)    # b^2 - 4*a*c
```

```@example poly
discriminant(x^3 - 3x + 2, x)        # repeated root at x = 1 → 0
```

```@example poly
discriminant(x^3 - 6x^2 + 11x - 6, x)  # roots 1, 2, 3 → 4
```

```@docs
Symbolics.discriminant
```

## Square-free decomposition

[`Symbolics.sqrfree`](@ref) decomposes a univariate polynomial into a product
of pairwise-coprime square-free factors with their multiplicities, via Yun's
algorithm. Given `f = c * ∏ p_i^{e_i}`, it returns the leading / content
scalar `c` and the list of pairs `(p_i, e_i)`.

```@example poly
u, fs = sqrfree((x - 1)^2 * (x - 2)^3, x)
```

The reconstruction `u · ∏ p_i^{e_i}` equals the input:

```@example poly
Symbolics.simplify(Symbolics.expand(u * prod(f^e for (f, e) in fs) - (x - 1)^2 * (x - 2)^3))
```

The leading content is extracted into `u`:

```@example poly
sqrfree(2x^2 - 2, x)
```

```@docs
Symbolics.sqrfree
```

## Factorisation over ℚ

[`Symbolics.factor`](@ref) produces the full factorisation of a univariate
polynomial over the rationals. It pulls out the integer content into the
unit and returns the remaining irreducible-over-ℚ factors with their
multiplicities.

`Symbolics.factor` is **not** exported, to avoid shadowing `Base.factor` and
`Primes.factor`; call it as `Symbolics.factor(...)`.

```@example poly
u, fs = Symbolics.factor(x^4 - 1, x)   # (x - 1)(x + 1)(x^2 + 1)
```

```@example poly
u, fs = Symbolics.factor(x^3 - 6x^2 + 11x - 6, x)   # roots 1, 2, 3
```

The leading content is returned in `u`:

```@example poly
u, fs = Symbolics.factor(2x^2 - 2, x)   # 2 * (x - 1)(x + 1)
```

Irreducible polynomials (over ℚ) come back as a single factor of
multiplicity 1:

```@example poly
Symbolics.factor(x^2 + 1, x)
```

### When the pure-Julia core isn't enough

The core pipeline is: pull out content via [`sqrfree`](@ref
Symbolics.sqrfree), then extract rational roots (yielding linear factors),
then classify any remaining residual. For a residual of degree at most 3
with no rational root, irreducibility is provable in pure Julia. For a
residual of degree 4 or more, the factorisation is best-effort; loading
`Nemo` delegates such residuals to `Nemo.factor` via the
`SymbolicsNemoExt` extension.

```@docs
Symbolics.factor
```
