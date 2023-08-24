using Test
using SymPy
using Symbolics

PythonCall.pyconvert(Symbolics.Num, sp.sympify("t"))

PythonCall.pyconvert(Symbolics.Num, sp.sympify("t**t"))

PythonCall.pyconvert(Symbolics.Num, sp.sympify("t + z + d"))

PythonCall.pyconvert(Symbolics.Num, sp.sympify("t*z*9"))

PythonCall.pyconvert(Symbolics.Num, sp.sympify("5*t*z + 3*d + h/(b*5)"))

PythonCall.pyconvert(Symbolics.Num, sp.sympify("t * n/z * t**4 * h**z + l*h - j"))

PythonCall.pyconvert(Symbolics.Equation, sp.sympify("Eq(2,5, evaluate = False)"))

PythonCall.pyconvert(Symbolics.Equation, sp.sympify("Eq(t*x + 5**x + 20/t, 90/t + t**4 - t*z)"))

PythonCall.pyconvert(Symbolics.Num, sp.sympify("Function('f')(x)"))

PythonCall.pyconvert(Symbolics.Num, sp.sympify("f(x,y)"))

PythonCall.pyconvert(Symbolics.Equation, sp.sympify("Eq(f(x), 2*x +1)"))

PythonCall.pyconvert(Symbolics.Num,sp.sympify("Derivative(f(x),x)"))

PythonCall.pyconvert(Symbolics.Num,sp.sympify("Eq(Derivative(f(x),x), x"))

fol = sp.sympify("Eq(Derivative(x(t),t), (1 - x(t))/3.0)")
convertedsystem = PythonCall.pyconvert(Symbolics.Equation,fol)

@named odesys = ODESystem(convertedsystem,t)

sys = structural_simplify(odesys)

prob = ODEProblem(sys, [0.0],(0.0,10.0))

using DifferentialEquations: solve
sol = solve(prob)

using Plots

plot(sol)

