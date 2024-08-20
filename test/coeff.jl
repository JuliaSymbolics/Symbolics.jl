using Symbolics, Test
import Symbolics: coeff

@variables x y z a b

@test isequal(coeff(2), 2)
@test isequal(coeff(2, x), 0)

@test isequal(coeff(2a, x), 0)
@test isequal(coeff(a*x, x), a)
@test isequal(coeff(2x*a, x), 2a)
# Symbolic powers:
@test isequal(coeff(a*x^b, x^b), a)
@test isequal(coeff(a*x^(b+1), x^(b+1)), a)
# Irrational powers:
@test isequal(coeff(a*x^sqrt(2), x^sqrt(2)), a)

@test isequal(coeff(a + x, x), 1)
@test isequal(coeff(2(a + x), x), 2)

e = 4 + x + 3x^2 + 2x^4 + a*x^2 + b
@test isequal(coeff(e), 4)
@test isequal(coeff(e, x^1), 1)
@test isequal(coeff(e, x^2), 3 + a)
@test isequal(coeff(e, x^3), 0)
@test isequal(coeff(e, x^4), 2)

e = x*y^2 + 2x + y^3*x^3
@test isequal(coeff(e, x), 2 + y^2)
@test isequal(coeff(e, x^3), y^3)
@test isequal(coeff(e, y^2), x)
@test isequal(coeff(e, y^3), x^3)

@test isequal(coeff(x^2 + y^2 + z^2, sin(z)), 0)
@test isequal(coeff(sin(z), z), 0)
@test isequal(coeff(sin(z)), 0)

# issue #236
@test isequal(coeff(3x + 2y, x), 3)
@test isequal(coeff(x*y, x), y)
@test isequal(coeff(x^2 + y, x^2), 1)

# expand - simplify needed
@test isequal(coeff(expand(((x + 1)^4 + x)^3), x^2), 93)
@test isequal(coeff(simplify((x^2 - 1) / (x - 1)), x), 1)
@test isequal(coeff(expand((x^(1//2) + y^0.5)^2), x), 1)


# issue #1098
@test isequal(coeff(x^2 + 1, x^0), 1)
@test isequal(coeff(e, x^0), 0)
@test isequal(coeff(a*x + 3, x^0), 3)

@test isequal(coeff(x / 5, x), 1//5)
@test isequal(coeff(x / y, x), 1/y)
@test isequal(coeff(x * 5y / (1 + y + z) , x), 5y / (1 + y + z))
