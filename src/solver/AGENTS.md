# Solver maintenance

Fraction clearing can introduce floating coefficients after the initial polynomial
filter. Normalize both cleared numerators and stored denominators before passing
them to Nemo or Groebner. Groebner's extension accepts systems as `Vector{Num}`.

Rational-solver regressions should cover univariate and multivariate inputs,
including single equations, wrapped and unwrapped expressions, floating
coefficients, and poles at irrational roots. Check hand-derived solutions and
substitution into the original equations. Run `test/solver.jl` in full and the
`GROUP=Core` test group after changing shared preprocessing.
