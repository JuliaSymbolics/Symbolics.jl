using Symbolics, Test

ex = [:(y ~ x)
      :(y ~ -2x + 3 / z)
      :(z ~ 2)]
eqs = parse_expr_to_symbolic.(ex, (@__MODULE__,))

@variables x y z
ex = [y ~ x
      y ~ -2x + 3 / z
      z ~ 2]
@test all(isequal.(eqs,ex))

ex = [:(b(t) ~ a(t))
      :(b(t) ~ -2a(t) + 3 / c(t))
      :(c(t) ~ 2)]
eqs = parse_expr_to_symbolic.(ex, (@__MODULE__,))
@variables t a(t) b(t) c(t)
ex = [b ~ a
      b ~ -2a + 3 / c
      c ~ 2]
@test_broken all(isequal.(eqs,ex))

# Unlike above tests variables need to be defined ahead of time
# To avoid BoundsError
@variables m[1:3]
ex = [:(m[2] ~ m[1])
      :(m[2] ~ -2m[1] + 3 / m[3])
      :(m[3] ~ 2)]
eqs = parse_expr_to_symbolic.(ex, (@__MODULE__,))
ex = [m[2] ~ m[1]
      m[2] ~ -2m[1] + 3 / m[3]
      m[3] ~ 2]
@test all(isequal.(eqs,ex))