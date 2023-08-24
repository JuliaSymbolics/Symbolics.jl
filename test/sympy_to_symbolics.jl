using Test
using PythonCall
using Symbolics

CondaPkg.add("sympy")

sp = pyimport("sympy")

@test Symbolics.tosymbol(sympy_to_symbolics(sp.sympify("t"))) == Symbol("t")

@test Symbolics.tosymbol(sympy_to_symbolics(sp.sympify("t**t"))) == Symbol("^(t, t)")

@test Symbolics.tosymbol(sympy_to_symbolics(sp.sympify("t + z + d"))) == Symbol("+(d, t, z)")

@test Symbolics.tosymbol(sympy_to_symbolics(sp.sympify("t*z*9"))) == Symbol("*(9, t, z)")

@test Symbolics.tosymbol(sympy_to_symbolics(sp.sympify("5*t*z + 3*d + h/(b*5)"))) == Symbol("+(3d, ((1//5)*h) / b, 5t*z)")

@test Symbolics.tosymbol(sympy_to_symbolics(sp.sympify("t * n/z * t**4 * h**z + l*h - j"))) == Symbol("+(-j, h*l, (n*(h^z)*(t^5)) / z)")

@test Symbolics.tosymbol(sympy_to_symbolics(sp.sympify("Eq(2,5, evaluate = False)")))

@test Symbol(sympy_to_symbolics(sp.sympify("Eq(t*x + 5**x + 20/t, 90/t + t**4 - t*z)"))) == Symbol("t^4 + 90 / t - t*z ~ 20 / t + t*x + 5^x")

@test Symbolics.tosymbol(sympy_to_symbolics(sp.sympify("Function('f')(x)"))) == Symbol("f(x)")

@test Symbolics.tosymbol(sympy_to_symbolics(sp.sympify("f(x,y)"))) == Symbol("f(x, y)")

@test Symbol((sympy_to_symbolics(sp.sympify("Eq(f(x), 2*x +1)")))) == Symbol("1 + 2x ~ f(x)")

@test Symbolics.tosymbol(sympy_to_symbolics(sp.sympify("Derivative(f(x),x)"))) == Symbol("fˍx(x)")

@test Symbol(sympy_to_symbolics(sp.sympify("Eq(Derivative(f(x),x), x)"))) == Symbol("x ~ Differential(x)(f(x))")
