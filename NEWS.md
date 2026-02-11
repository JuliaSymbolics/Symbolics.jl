## 7.0.0

### Breaking changes

- `substitute` no longer recurses into `Differential` arguments. This is due to
  SymbolicUtils.jl v4's `default_substitute_filter`, which treats `Operator` subclasses
  (including `Differential`) as substitution boundaries. Use the new `substitute_in_deriv`
  or `substitute_in_deriv_and_depvar` functions to substitute inside `Differential`
  expressions. See the [Derivatives documentation](https://docs.sciml.ai/Symbolics/stable/manual/derivatives/)
  for details.

## 4.0.0
- Symbolics.jl now supports the latest symbolic computing architecture backed by Metatheory.jl v1.2
  and SymbolicUtils.jl v0.18 for generic term rewriting.
- Support for automatic code optimization through Metatheory.jl EGraphs and SymbolicUtils's `optimize` function.

## 3.3.0

- adds `simplify_fractions` which turns an expression into a single fraction
  and simplifies by dividing the numerator and denominator factors by
  appropriate GCDs
- Use new `fraction_iszero` and `fraction_isone` functions from SymbolicUtils
  to implement `iszero` and `isone` respectively.
- `x / x` etc. are no more simplified on construction, call
  `simplify_fractions` to simplify them.
