using Symbolics, Test

ex = [:(y ~ x)
      :(y ~ -2x + 3 / z)
      :(z ~ 2)]
eqs = parse_expr_to_symbolic.(ex, (Main,))

@variables x y z
ex = [y ~ x
      y ~ -2x + 3 / z
      z ~ 2]
@test all(isequal.(eqs,ex))

ex = [:(b(t) ~ a(t))
      :(b(t) ~ -2a(t) + 3 / c(t))
      :(c(t) ~ 2)]
eqs = parse_expr_to_symbolic.(ex, (Main,))
@variables t a(t) b(t) c(t)
ex = [b ~ a
      b ~ -2a + 3 / c
      c ~ 2]
@test_broken all(isequal.(eqs,ex))