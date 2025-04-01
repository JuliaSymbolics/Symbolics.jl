using Symbolics, Test
using DomainSets
@variables x y

I1 = Integral(x in ClosedInterval(1, 5))
I2 = Integral(x in ClosedInterval(1, 5))
@test I1 == I2

@variables v(..) u(..) x y a b

# test constant integrand
I = Integral(x in ClosedInterval(a, b))
@test isequal(I(0), 0)
@test isequal(I(2), 2*(b -a))
@test isequal(I(2.1), 2.1*(b - a))
@test isequal(I(pi), pi*(b - a))
@test isequal(I(1//2), 1//2 * (b - a))

# test complex integrand
@test I(2im) isa Complex{Num}
@test isequal(I(2im), 2im * (b - a))
@test isequal(I(1 + 2.1im), (1 + 2.1im)*(b - a))
@test I(x + imx) isa Complex{Num}

D = Differential(x)
Dxx = Differential(x)^2
I = Integral(y in ClosedInterval(1, 5))
eq = D(I(u(x,y))) ~ 0
eq_test = I(D(u(x,y)))
@test isequal(expand_derivatives(eq.lhs), Symbolics.value(eq_test))

eq = Dxx((I(u(x,y)))) + I(D(u(x,y))) ~ 0
eq_test = I(D(u(x,y))) + I(D(D(u(x,y))))
@test isequal(expand_derivatives(eq.lhs), Symbolics.value(eq_test))

I = Integral(y in ClosedInterval(1, v(x)))
eq = D((I(u(x,y)))) ~ 0
eq_test = I(D(u(x,y))) + D(v(x))*u(x, v(x))
@test isequal(expand_derivatives(eq.lhs), Symbolics.value(eq_test))

I = Integral(y in ClosedInterval(1, v(x)))
eq = Dxx((I(u(x,y)))) ~ 0
eq_test = D(I(D(u(x,y))) + D(v(x))*u(x, v(x)))
eq_test_ = I(D(D(u(x,y)))) + D(D(v(x)))*u(x, v(x)) + 2D(u(x,v(x)))*D(v(x))
@test isequal(expand_derivatives(eq.lhs), Symbolics.value(eq_test_))
@test isequal(expand_derivatives(eq.lhs), expand_derivatives(eq_test))

eq = D((I(u(x,y)^2))) ~ 0
eq_test = I(2D(u(x,y))*u(x,y)) + D(v(x))*u(x, v(x))^2
@test isequal(expand_derivatives(eq.lhs), Symbolics.value(eq_test))

# case where limits of integral contain the variable to derive
# against
I = Integral(y in ClosedInterval(1, 2x))
@test isequal(expand_derivatives(D(I(1))), 2)

# same but case where limit of integral is not a call
I = Integral(y in ClosedInterval(1, x))
@test isequal(expand_derivatives(D(I(1))), 1)

# test shadowing by integration variable 
# the result of a definite integral over x does not depend on x
# anymore unless it appears again in the limits
I = Integral(x in ClosedInterval(1, 2))
@test isequal(expand_derivatives(D(I(u(x)))), 0)
