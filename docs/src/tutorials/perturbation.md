# [Mixed Symbolic-Numeric Perturbation Theory](@id perturb_alg)

[**Symbolics.jl**](https://github.com/JuliaSymbolics/Symbolics.jl) is a fast and modern Computer Algebra System (CAS) written in Julia. It is an integral part of the [SciML](https://sciml.ai/) ecosystem of differential equation solvers and scientific machine learning packages. While **Symbolics.jl** is primarily designed for modern scientific computing (e.g. automatic differentiation and machine learning), it is also a powerful CAS that can be used for *classic* scientific computing. One such application is using *perturbation theory* to solve algebraic and differential equations.

Perturbation methods are a collection of techniques to solve hard problems that generally don't have a closed solution, but depend on a tunable parameter and have closed or easy solutions for some values of this parameter. The main idea is to assume a solution that is a power series in the tunable parameter (say $ϵ$), such that $ϵ = 0$ corresponds to an easy solution, and then solve iteratively for higher-order corrections.

The hallmark of perturbation theory is the generation of long and convoluted intermediate equations by this process. These are subjected to algorithmic and mechanical manipulations, which makes perturbation methods an excellent fit for a CAS. In fact, CAS software have been used for perturbation calculations since the early 1970s.

This tutorial shows how to mix symbolic manipulations and numerical methods to solve algebraic equations with perturbation theory. [Another tutorial applies it to differential equations](https://docs.sciml.ai/ModelingToolkit/stable/examples/perturbation/).

## Solving the quintic equation

The “hello world!” analog of perturbation problems is to find a real solution $x$ to the quintic (fifth-order) equation
```@example perturb
using Symbolics
@variables x
quintic = x^5 + x ~ 1
```
According to Abel's theorem, a general quintic equation does not have a closed form solution. But we can easily solve it numerically using Newton's method (here implemented for simplicity, and not performance):
```@example perturb
function solve_newton(eq, x, x₀; abstol=1e-8, maxiters=50)
    # symbolic expressions for f(x) and f′(x)
    f = eq.lhs - eq.rhs # want to find root of f(x)
    f′ = Symbolics.derivative(f, x)

    xₙ = x₀ # numerical value of the initial guess
    for i = 1:maxiters
        # calculate new guess by numerically evaluating symbolic expression at previous guess
        xₙ₊₁ = substitute(x - f / f′, x => xₙ)
        if abs(xₙ₊₁ - xₙ) < abstol
            return xₙ₊₁ # converged
        else
            xₙ = xₙ₊₁
        end
    end
    error("Newton's method failed to converge")
end

x_newton = solve_newton(quintic, x, 1.0)
println("Newton's method solution: x = ", x_newton)
```

To solve the same problem with perturbation theory, we must introduce an expansion variable $ϵ$ in the equation:
```@example perturb
@variables ϵ # expansion variable
quintic = x^5 + ϵ*x ~ 1
```
If $ϵ = 1$, we get our original problem. With $ϵ = 0$, the problem transforms to the easy quintic equation $x^5 = 1$ with the trivial real solution $x = 1$ (and four complex solutions which we ignore). Next, expand $x$ as a power series in $ϵ$:
```@example perturb
x_taylor = series(x, ϵ, 0:7) # expand x in a power series in ϵ
x_coeffs = taylor_coeff(x_taylor, ϵ) # TODO: get coefficients at series creation
x_taylor
```
Then insert this into the quintic equation and expand it, too, to the same order:
```@example perturb
quintic_taylor = substitute(quintic, x => x_taylor)
quintic_taylor = taylor(quintic_taylor, ϵ, 0:7)
```
This messy equation must hold for each power of $ϵ$, so we can separate it into one nicer equation per order:
```@example perturb
quintic_eqs = taylor_coeff(quintic_taylor, ϵ)
quintic_eqs[1:5] # for readability, show only 5 shortest equations
```
These equations show three important features of perturbation theory:
1. The $0$-th order equation is *trivial* in $x_0$: here $x_0^5 = 1$ has the trivial real solution $x_0 = 1$.
2. The $n$-th order equation is *linear* in $x_n$ (except the trivial $0$-th order equation).
3. The equations are *triangular* in $x_n$: the $n$-th order equation can be solved for $x_n$ given only $x_m$ for $m<n$.
This structure is what makes the perturbation theory so attractive: we can start with the trivial solution $x_0 = 1$, then linearly solve for $x_n$ order-by-order and substitute in the solutions of $x_m$ for $m<n$ obtained so far.

Here is a simple function that uses this *cascading* strategy to solve such a set of equations `eqs` for the variables `xs`, given a solution `x₀` of the first equation `eqs[begin]`:
```@example perturb
function solve_cascade(eqs, xs, x₀, ϵ)
    sol = Dict(xs[begin] => x₀) # store solutions in a map

    # verify that x₀ is a solution of the first equation
    eq0 = substitute(eqs[1], sol)
    isequal(eq0.lhs, eq0.rhs) || error("$sol does not solve $(eqs[1])")

    # solve remaining equations sequentially
    for i in 2:lastindex(eqs)
        eq = substitute(eqs[i], sol) # insert previous solutions
        x = Symbolics.symbolic_linear_solve(eq, xs[i]) # solve current equation
        sol = merge(sol, Dict(xs[i] => x)) # store solution
    end

    return sol
end
```
Let us solve our order-separated quintics for the coefficients, and substitute them into the full series for $x$:
```@example perturb
x_coeffs_sol = solve_cascade(quintic_eqs, x_coeffs, 1, ϵ)
x_pert = substitute(x_taylor, x_coeffs_sol)
```
The $n$-th order solution of our original quintic equation is the sum up to the $ϵ^n$-th order term, evaluated at $ϵ=1$:
```@example perturb
for n in 0:7
    x_pert_sol = substitute(taylor(x_pert, ϵ, 0:n), ϵ => 1)
    println("$n-th order solution: x = $x_pert_sol = $(x_pert_sol * 1.0)")
end
```
This is close to the solution from Newton's method!

## Solving Kepler's Equation

Historically, perturbation methods were first invented to calculate orbits of the Moon and the planets. In homage to this history, our second example is to solve [Kepler's equation](https://en.wikipedia.org/wiki/Kepler's_equation), which is central to solving two-body Keplerian orbits:
```@example perturb
@variables e E M
kepler = E - e * sin(E) ~ M
```
We want to solve for the *eccentric anomaly* $E$ given the *eccentricity* $e$ and *mean anomaly* $M$.
This is also easy with Newton's method. With Earth's eccentricity $e = 0.01671$ and $M = \pi/2$:
```@example perturb
vals_earth = Dict(e => 0.01671, M => π/2)
E_newton = solve_newton(substitute(kepler, vals_earth), E, π/2)
println("Newton's method solution: E = ", E_newton)
```

Next, let us solve the same problem with the perturbative method. It is most common to expand $E$ as a series in $M$. Repeating the procedure from the quintic example, we get these equations:
```@example perturb
E_taylor = series(E, M, 0:5)
E_coeffs = taylor_coeff(E_taylor, M) # TODO: get coefficients at series creation
kepler_eqs = taylor_coeff(substitute(kepler, E => E_taylor), M, 0:5)
kepler_eqs[1:4] # for readability
```
The trivial $0$-th order solution (when $M=0$) is $E_0=0$. This gives this full perturbative solution:
```@example perturb
E_coeffs_sol = solve_cascade(kepler_eqs, E_coeffs, 0, M)
E_pert = substitute(E_taylor, E_coeffs_sol)
```

Numerically, the result again converges to that from Newton's method:
```@example perturb
for n in 0:5
    println("$n-th order solution: E = ", substitute(taylor(E_pert, M, 0:n), vals_earth))
end
```
But unlike Newtons method, this example shows how perturbation theory also gives us the powerful *symbolic* series solution for $E$ (*before* numbers for $e$ and $M$ are inserted). Our series matches [this result from Wikipedia](https://en.wikipedia.org/wiki/Kepler%27s_equation#Inverse_Kepler_equation):
```@example perturb
E_wiki = 1/(1-e)*M - e/(1-e)^4*M^3/factorial(3) + (9e^2+e)/(1-e)^7*M^5/factorial(5)
```

Alternatively, we can expand $E$ in $e$ instead of $M$, giving the solution (the trivial solution when $e = 0$ is $E_0=M$):
```@example perturb
E_taylor′ = series(E, e, 0:5)
E_coeffs′ = taylor_coeff(E_taylor′, e) # TODO: get at creation
kepler_eqs′ = taylor_coeff(substitute(kepler, E => E_taylor′), e, 0:5)
E_coeffs_sol′ = solve_cascade(kepler_eqs′, E_coeffs′, M, e)
E_pert′ = substitute(E_taylor′, E_coeffs_sol′)
```
This looks very different from our first series `E_pert`. If they are the same, we should get $0$ if we subtract and expand both as multivariate Taylor series in $(e,M)$. Indeed:
```@example perturb
@assert taylor(taylor(E_pert′ - E_pert, e, 0:4), M, 0:4) == 0 # use this as a test # hide
taylor(taylor(E_pert′ - E_pert, e, 0:4), M, 0:4)
```
