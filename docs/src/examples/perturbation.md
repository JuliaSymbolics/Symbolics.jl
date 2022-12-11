# [Mixed Symbolic-Numeric Perturbation Theory](@id perturb_alg)

## Background

[**Symbolics.jl**](https://github.com/JuliaSymbolics/Symbolics.jl) is a fast and modern Computer Algebra System (CAS) written in the Julia Programming Language. It is an integral part of the [SciML](https://sciml.ai/) ecosystem of differential equation solvers and scientific machine learning packages. While **Symbolics.jl** is primarily designed for modern scientific computing (e.g., auto-differentiation, machine learning), it is a powerful CAS and can also be useful for *classic* scientific computing. One such application is using the *perturbation* theory to solve algebraic and differential equations.

Perturbation methods are a collection of techniques to solve intractable problems that generally don't have a closed solution but depend on a tunable parameter and have closed or easy solutions for some values of the parameter. The main idea is to assume a solution as a power series in the tunable parameter (say $ϵ$), such that $ϵ = 0$ corresponds to an easy solution.

We will discuss the general steps of the perturbation methods to solve algebraic (this tutorial) and [differential equations](@ref perturb_diff)

The hallmark of the perturbation method is the generation of long and convoluted intermediate equations, which are subjected to algorithmic and mechanical manipulations. Therefore, these problems are well suited for CAS. In fact, CAS softwares have been used to help with the perturbation calculations since the early 1970s.

In this tutorial our goal is to show how to use a mix of symbolic manipulations (**Symbolics.jl**) and numerical methods to solve simple perturbation problems.

## Solving the Quintic

We start with the "hello world!" analog of the perturbation problems, solving the quintic (fifth-order) equations. We want to find a real valued $x$ such that $x^5 + x = 1$. According to the Abel's theorem, a general quintic equation does not have a closed form solution. Of course, we can easily solve this equation numerically; for example, by using the Newton's method. We use the following implementation of the Newton's method:

```@example perturb
using Symbolics, SymbolicUtils

function solve_newton(f, x, x₀; abstol=1e-8, maxiter=50)
    xₙ = Float64(x₀)
    fₙ₊₁ = x - f / Symbolics.derivative(f, x)

    for i = 1:maxiter
        xₙ₊₁ = substitute(fₙ₊₁, Dict(x => xₙ))
        if abs(xₙ₊₁ - xₙ) < abstol
            return xₙ₊₁
        else
            xₙ = xₙ₊₁
        end
    end
    return xₙ₊₁
end
```

In this code, `Symbolics.derivative(eq, x)` does exactly what it names implies: it calculates the symbolic derivative of `eq` (a **Symbolics.jl** expression) with respect to `x` (a **Symbolics.jl** variable). We use `Symbolics.substitute(eq, D)` to evaluate the update formula by substituting variables or sub-expressions (defined in a dictionary `D`) in `eq`. It should be noted that `substitute` is the workhorse of our code and will be used multiple times in the rest of these tutorials. `solve_newton` is written with simplicity and clarity, and not performance, in mind but suffices for our purpose.

Let's go back to our quintic. We can define a Symbolics variable as `@variables x` and then solve the equation `solve_newton(x^5 + x - 1, x, 1.0)` (here, `x₀ = 0` is our first guess). The answer is 0.7549. Now, let's see how we can solve the same problem using the perturbation methods.

We introduce a tuning parameter $\epsilon$ into our equation: $x^5 + \epsilon x = 1$. If $\epsilon = 1$, we get our original problem. For $\epsilon = 0$, the problem transforms to an easy one: $x^5 = 1$ which has an exact real solution $x = 1$ (and four complex solutions which we ignore here). We expand $x$ as a power series on $\epsilon$:

```math
  x(\epsilon) = a_0 + a_1 \epsilon + a_2 \epsilon^2 + O(\epsilon^3)
```

$a_0$ is the solution of the easy equation, therefore $a_0 = 1$. Substituting into the original problem,

```math
  (a_0 + a_1 \epsilon + a_2 \epsilon^2)^5 + \epsilon (a_0 + a_1 \epsilon + a_2 \epsilon^2) - 1 = 0
```

Expanding the equations, we get

```math
  \epsilon (1 + 5 a_1) + \epsilon^2 (a_1 + 5 a_2 + 10 a1_2) + 𝑂(\epsilon^3) = 0
```

This equation should hold for each power of $\epsilon$. Therefore,

```math
  1 + 5 a_1 = 0
```

and

```math
  a_1 + 5 a_2 + 10 a_1^2 = 0
```

This system of equations does not initially seem to be linear because of the presence of terms like $10 a_1^2$, but upon closer inspection is found to be in fact linear (this is a feature of the perturbation methods). In addition, the system is in a triangular form, meaning the first equation depends only on $a_1$, the second one on $a_1$ and $a_2$, such that we can replace the result of $a_1$ from the first one into the second equation and remove the non-linear term. We solve the first equation to get $a_1 = -\frac{1}{5}$. Substituting in the second one and solve for $a_2$:

```math
  a_2 = \frac{(-\frac{1}{5} + 10(-(\frac{1}{5})²)}{5}  = -\frac{1}{25}
```

Finally,

```math
  x(\epsilon) = 1 - \frac{\epsilon}{5} - \frac{\epsilon^2}{25} + O(\epsilon^3)
```

Solving the original problem, $x(1) = 0.76$, compared to 0.7548 calculated numerically. We can improve the accuracy by including more terms in the expansion of $x$. However, the calculations, while straightforward, become messy and intractable to do manually very quickly. This is why a CAS is very helpful to solve perturbation problems.

Now, let's see how we can do these calculations in Julia. Let $n$ be the order of the expansion. We start by defining the symbolic variables:

```@example perturb
n = 2
@variables ϵ a[1:n]
```

Then, we define

```@example perturb
x = 1 + a[1]*ϵ + a[2]*ϵ^2
```

The next step is to substitute `x` in the problem equation

```@example perturb
  eq = x^5 + ϵ*x - 1
```

The expanded form of `eq` is

```@example perturb
expand(eq)
```

We need a way to get the coefficients of different powers of `ϵ`. Function `collect_powers(eq, x, ns)` returns the powers of variable `x` in expression `eq`. Argument `ns` is the range of the powers.

```@example perturb
function collect_powers(eq, x, ns; max_power=100)
    eq = substitute(expand(eq), Dict(x^j => 0 for j=last(ns)+1:max_power))

    eqs = []
    for i in ns
        powers = Dict(x^j => (i==j ? 1 : 0) for j=1:last(ns))
        push!(eqs, substitute(eq, powers))
    end
    eqs
end
```

To return the coefficients of $ϵ$ and $ϵ^2$ in `eq`, we can write

```@example perturb
eqs = collect_powers(eq, ϵ, 1:2)
```

A few words on how `collect_powers` works, It uses `substitute` to find the coefficient of a given power of `x` by passing a `Dict` with all powers of `x` set to 0, except the target power which is set to 1. For example, the following expression returns the coefficient of `ϵ^2` in `eq`,

```@example perturb
substitute(expand(eq), Dict(
  ϵ => 0,
  ϵ^2 => 1,
  ϵ^3 => 0,
  ϵ^4 => 0,
  ϵ^5 => 0,
  ϵ^6 => 0,
  ϵ^7 => 0,
  ϵ^8 => 0)
)
```

Back to our problem. Having the coefficients of the powers of `ϵ`, we can set each equation in `eqs` to 0 (remember, we rearrange the problem such that `eq` is 0) and solve the system of linear equations to find the numerical values of the coefficients. **Symbolics.jl** has a function `Symbolics.solve_for` that can solve systems of linear equations. However, the presence of higher order terms in `eqs` prevents `Symbolics.solve_for(eqs .~ 0, a)` from workings properly. Instead, we can exploit the fact that our system is in a triangular form and start by solving `eqs[1]` for `a₁` and then substitute this in `eqs[2]` and solve for `a₂` (as continue the same process for higher order terms).  This *cascading* process is done by function `solve_coef(eqs, ps)`:

```@example perturb
function solve_coef(eqs, ps)
    vals = Dict()

    for i = 1:length(ps)
        eq = substitute(eqs[i], vals)
        vals[ps[i]] = Symbolics.solve_for(eq ~ 0, ps[i])
    end
    vals
end
```

Here, `eqs` is an array of expressions (assumed to be equal to 0) and `ps` is an array of variables. The result is a dictionary of *variable* => *value* pairs. We apply `solve_coef` to `eqs` to get the numerical values of the parameters:

```@example perturb
solve_coef(eqs, a)
```

Finally, we substitute back the values of `a` in the definition of `x` as a function of `𝜀`. Note that `𝜀` is a number (usually Float64), whereas `ϵ` is a symbolic variable.

```@example perturb
X = 𝜀 -> 1 + a[1]*𝜀 + a[2]*𝜀^2
```

Therefore, the solution to our original problem becomes `X(1)`, which is equal to 0.76. We can use larger values of `n` to improve the accuracy of estimations.

| n | x              |
|---|----------------|
|1  |0.8 |
|2  |0.76|
|3  |0.752|
|4  |0.752|
|5  |0.7533|
|6  |0.7543|
|7  |0.7548|
|8  |0.7550|

Remember the numerical value is 0.7549. The two functions `collect_powers` and `solve_coef(eqs, a)` are used in all the examples in this and the next tutorial.

## Solving the Kepler's Equation

Historically, the perturbation methods were first invented to solve orbital calculations of the Moon and the planets. In homage to this history, our second example has a celestial theme. Our goal is solve the Kepler's equation:

```math
  E - e\sin(E) = M
```

where $e$ is the *eccentricity* of the elliptical orbit, $M$ is the *mean anomaly*, and $E$ (unknown) is the *eccentric anomaly* (the angle between the position of a planet in an elliptical orbit and the point of periapsis). This equation is central to solving two-body Keplerian orbits.

Similar to the first example, it is easy to solve this problem using the Newton's method. For example, let $e = 0.01671$ (the eccentricity of the Earth) and $M = \pi/2$. We have `solve_newton(x - e*sin(x) - M, x, M)` equals to 1.5875 (compared to π/2 = 1.5708). Now, we try to solve the same problem using the perturbation techniques (see function `test_kepler`).

For $e = 0$, we get $E = M$. Therefore, we can use $e$ as our perturbation parameter. For consistency with other problems, we also rename $e$ to $\epsilon$ and $E$ to $x$.

From here on, we use the helper function `def_taylor` to define Taylor's series by calling it as `x = def_taylor(ϵ, a, 1)`, where the arguments are, respectively, the perturbation variable, an array of coefficients (starting from the coefficient of $\epsilon^1$), and an optional constant term.

```@example perturb
def_taylor(x, ps) = sum([a*x^i for (i,a) in enumerate(ps)])
def_taylor(x, ps, p₀) = p₀ + def_taylor(x, ps)
```

We start by defining the variables (assuming `n = 3`):

```@example perturb
n = 3
@variables ϵ M a[1:n]
x = def_taylor(ϵ, a, M)
```

We further simplify by substituting `sin` with its power series using the `expand_sin` helper function:

```@example perturb
expand_sin(x, n) = sum([(isodd(k) ? -1 : 1)*(-x)^(2k-1)/factorial(2k-1) for k=1:n])
```

To test,

```@example perturb
expand_sin(0.1, 10) ≈ sin(0.1)
```

The problem equation is

```@example perturb
eq = x - ϵ * expand_sin(x, n) - M
```

We follow the same process as the first example. We collect the coefficients of the powers of `ϵ`

```@example perturb
eqs = collect_powers(eq, ϵ, 1:n)
```

and then solve for `a`:

```@example perturb
vals = solve_coef(eqs, a)
```

Finally, we substitute `vals` back in `x`:

```@example perturb
x′ = substitute(x, vals)
X = (𝜀, 𝑀) -> substitute(x′, Dict(ϵ => 𝜀, M => 𝑀))
X(0.01671, π/2)
```

The result is 1.5876, compared to the numerical value of 1.5875. It is customary to order `X` based on the powers of `𝑀` instead of `𝜀`. We can calculate this series as `collect_powers(sol, M, 0:3)
`. The result (after manual cleanup) is

```
(1 + 𝜀 + 𝜀^2 + 𝜀^3)*𝑀
- (𝜀 + 4*𝜀^2 + 10*𝜀^3)*𝑀^3/6
+ (𝜀 + 16*𝜀^2 + 91*𝜀^3)*𝑀^5/120
```

Comparing the formula to the one for 𝐸 in the [Wikipedia article on the Kepler's equation](https://en.wikipedia.org/wiki/Kepler%27s_equation):

```math
  E = \frac{1}{1-\epsilon}M
    -\frac{\epsilon}{(1-\epsilon)^4} \frac{M^3}{3!} + \frac{(9\epsilon^2
    + \epsilon)}{(1-\epsilon)^7}\frac{M^5}{5!}\cdots
```

The first deviation is in the coefficient of $\epsilon^3 M^5$.
