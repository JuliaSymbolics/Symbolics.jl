# Partial Fraction Decomposition

Partial fraction decomposition is performed using the cover-up method. This involves "covering up" a factor in the denominator and substituting the root into the remaining expression. When the denominator can be completely factored into non-repeated linear factors, this produces the desired result. When there are repeated or irreducible quadratic factors, it produces terms with unknown coefficients in the numerator that is solved as a system of equations.

It is often used when solving integrals or performing an inverse Laplace transform (see [`inverse_laplace`](@ref)).

```docs
Symbolics.partial_frac_decomposition
```