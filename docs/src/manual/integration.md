# Symbolic Integrals

Symbolics.jl provides a the `Integral` operator for defining integrals. These can be used with
various functionalities in order to represent integro-differential equations and as a representation
for symbolic solving of integrals.

!!! note
    This area is currently under heavy development. More solvers will be available in the near future.

## Defining Symbolic Integrals

Note that integration domains are defined using 
[DomainSets.jl](https://github.com/JuliaApproximation/DomainSets.jl).

```@docs
Symbolics.Integral
```

## Solving Symbolic Integrals

Symbolics.jl currently has the following options for solving integrals symbolically:

| Option    | Description | Pros | Cons |
| -------- | ------- | ------- | ------- |
| SymbolicNumericIntegration.jl  | Uses numeric integrators with symbolic regression machine learning methods | Can solve some hard integrals very fast and easily | Can be unreliable in easy cases, will give floats (0.5) instead of rational values (1//2) |
| SymPy's Integrate | Uses SymPy (through SymPy.jl) to integrate the expression | Reasonably robust and tested solution | Extremely slow |

### SymbolicNumericIntegration.jl

SymbolicNumericIntegration.jl is a package for solving symbolic integrals using numerical methods. 
Specifically, it mixes numerical integration with machine learning symbolic regression in order to
discover the basis for the integrals solution, for which it then finds the coefficients and uses
differentiation to prove the correctness of the result. While this technique is unusual and new,
[see the paper for details](https://arxiv.org/abs/2201.12468v2), it has the advantage of working
very differently from the purely symbolic methods, and thus the problems which it finds hard can
be entirely different from ones which the purely symbolic methods find hard. That is, while purely
symbolic methods have difficulty depending on the complexity of the expression, this method has
difficulty depending on the uniqueness of the numerical solution of the expression, and there are many
complex expressions which have a very distinct solution behavior and thus can be easy to solve through
this method.

One major downside to this method is that, because it works numerically, there can be precision loss
in the coefficients. Thus the returned solution may have `0.5` instead of `1//2`. 

For using SymbolicNumericIntegration.jl, see [the SymbolicNumericIntegration.jl documentation](https://docs.sciml.ai/SymbolicNumericIntegration/stable/) for more details. A quick example is:

```julia
using Symbolics
using SymbolicNumericIntegration

@variables x a b

integrate(3x^3 + 2x - 5)
# (x^2 + (3 // 4) * (x^4) - (5 // 1) * x, 0, 0)

integrate(exp(a * x), x; symbolic = true)
# (exp(a * x) / a, 0, 0)

integrate(sin(a * x) * cos(b * x), x; symbolic = true, detailed = false)
# (-a * cos(a * x) * cos(b * x) - b * sin(a * x) * sin(b * x)) / (a^2 - (b^2))
```

### SymPy

The function `sympy_integrate` allows one to call 
[SymPy's integration functionality](https://docs.sympy.org/latest/modules/integrals/integrals.html).
While it's generally slow, if you leave your computer on for long enough it can solve some difficult
expressions.

```@docs
Symbolics.sympy_integrate
```
