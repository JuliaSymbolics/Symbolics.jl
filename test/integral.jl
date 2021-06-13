using Symbolics, Test
using DomainSets
@variables x y

I1 = Integral(x, ClosedInterval(1, 5))
I2 = Integral(x, ClosedInterval(1, 5))
@test I1 == I2

@variables v(..) u(..) x y
D = Differential(x)
Dxx = Differential(x)^2
I = Integral(y, ClosedInterval(1, 5))
eq = D((I(u(x,y)))) ~ 0
eq_test = I(D(u(x,y)))
@test isequal(expand_derivatives(eq.lhs), Symbolics.value(eq_test))

eq = Dxx((I(u(x,y)))) + I(D(u(x,y))) ~ 0
eq_test = I(D(u(x,y))) + I(D(D(u(x,y))))
@test isequal(expand_derivatives(eq.lhs), Symbolics.value(eq_test))

I = Integral(y, ClosedInterval(1, v(x)))
eq = D((I(u(x,y)))) ~ 0
eq_test = I(D(u(x,y))) + D(v(x))*u(x, v(x))
@test isequal(expand_derivatives(eq.lhs), Symbolics.value(eq_test)) == true

I = Integral(y, ClosedInterval(1, v(x)))
eq = Dxx((I(u(x,y)))) ~ 0
eq_test = D(I(D(u(x,y))) + D(v(x))*u(x, v(x)))
eq_test_ = I(D(D(u(x,y)))) + D(D(v(x)))*u(x, v(x)) + 2D(u(x,v(x)))*D(v(x))
@test isequal(expand_derivatives(eq.lhs), Symbolics.value(eq_test_)) == true
@test isequal(expand_derivatives(eq.lhs), expand_derivatives(eq_test)) == true
