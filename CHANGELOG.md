# Changelog

All notable changes to Symbolics.jl will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [7.20.0] - 2026-04-21

### Added

- `resultant(f, g, var; algorithm = :euclid | :sylvester)` for univariate
  polynomials with rational or symbolic coefficients. Both algorithms are
  exposed under a single entry point and produce mathematically equal results
  on every valid input.
- `discriminant(f, var)` for univariate polynomials, derived from the
  classical identity `disc(f) = (-1)^(n(n-1)/2) * Res(f, f') / lc(f)`.
  Returns `0` when `f` has a repeated root, non-zero otherwise. Conventions:
  `discriminant(a*x + b, x) == 1` for linear inputs; `discriminant(c, x) == 0`
  for a non-zero constant.
- `sqrfree(f, var)` for square-free decomposition of univariate polynomials
  via Yun's algorithm (characteristic 0). Returns `(unit, factors)` where
  `factors::Vector{Tuple{Num, Int}}` is a list of `(p_i, e_i)` pairs with
  pairwise-coprime square-free `p_i`; `unit * prod(p_i^e_i) == f`.
- `Symbolics.factor(f, var)` for univariate factorisation over ℚ. Pulls
  content into the unit, applies square-free decomposition, and extracts
  rational linear factors; degree-≤-3 residuals with no rational root are
  provably irreducible. Not exported to avoid shadowing `Base.factor` /
  `Primes.factor`. Higher-degree residuals delegate to `Nemo.factor` via
  the `SymbolicsNemoExt` extension when `Nemo` is loaded (fully wired:
  `factor(x^10 - 1, x)` returns the four cyclotomic factors over ℚ with
  Nemo loaded, and only the partial three-factor split without it).
  Floating-point coefficients are rejected with a clear error.

### Changed

- `SymbolicsNemoExt.nemo_crude_evaluate`: added explicit dispatches for
  `Nemo.QQFieldElem` and tightened the `Nemo.FracElem` method so that
  integer-valued rationals (`n//1`) collapse to plain `Int`/`BigInt`
  constants in returned expressions. Purely cosmetic — the existing
  `Symbolics.factor_use_nemo` (used by the solver) still produces
  mathematically identical output but without the noisy `1//1` literals.
