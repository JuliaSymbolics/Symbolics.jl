## 3.3.0

- adds `simplify_fractions` which turns an expression into a single fraction
  and simplifies by dividing the numerator and denominator factors by
  appropriate GCDs
- Use new `fraction_iszero` and `fraction_isone` functions from SymbolicUtils
  to implement `iszero` and `isone` respectively.
- `x / x` etc. are no more simplified on construction, call
  `simplify_fractions` to simplify them.

## 3.4.0

- adds `optimize` function which can optimize an expression to minimize CPU usage
  by re-arranging the expression. It does not guarrantee numerically identical results.
- allows you to use `@methodtheory` macro to define equational rules using Metatheory.jl
