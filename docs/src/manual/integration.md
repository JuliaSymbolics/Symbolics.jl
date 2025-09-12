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

| Package    | Method | Description | Pros | Cons |
| ---------- | ------ | ----------- | ---- | ---- |
| SymbolicIntegration.jl | Risch Algorithm | Implements theoretical Risch algorithm for elementary functions | Mathematically rigorous, returns exact symbolic results, handles rational functions and transcendental functions (exp, log, trig) | Still in early development, unstable, limited to specific function classes, no support for algebraic functions |
| SymbolicIntegration.jl | Rule-Based Method | Applies over 3,400 integration rules systematically | Handles non-integer powers, supports multiple symbolic variables, more flexible than pure Risch | Rule-dependent success, no theoretical completeness guarantee, limited by predefined rule set |
| SymbolicNumericIntegration.jl  | Symbolic-Numeric Hybrid | Uses numeric integrators with symbolic regression machine learning methods | Can solve some hard integrals very fast and easily, works differently from pure symbolic methods | Can be unreliable in easy cases, will give floats (0.5) instead of rational values (1//2), precision loss |
| SymPy's Integrate | Hybrid Methods | Uses SymPy (through SymPy.jl) to integrate the expression | Reasonably robust and tested solution, mature implementation | Extremely slow, Python overhead |

### SymbolicIntegration.jl

SymbolicIntegration.jl provides pure symbolic integration using multiple algorithmic approaches.
Unlike numerical methods, it returns exact symbolic results with rational coefficients (e.g., `1//2`
instead of `0.5`). The package implements a flexible framework for symbolic integration that uses Julia's
multiple dispatch to select appropriate algorithms for different types of integrands.

The package offers two main integration methods:

#### Risch Algorithm Method

The `RischMethod` implements the theoretical Risch algorithm for integrating:
- Polynomial functions
- Rational functions  
- Exponential functions
- Logarithmic functions
- Trigonometric functions
- Combinations of the above

This method is mathematically rigorous and provides guaranteed results when an elementary antiderivative exists.
However, it is still under development and may not handle all cases robustly.

#### Rule-Based Method

The rule-based approach systematically applies over 3,400 integration rules to solve integrals:
- Handles non-integer powers and algebraic functions
- Supports more complex integrations than pure Risch
- Can handle multiple symbolic variables
- More flexible but depends on the predefined rule set

The package automatically tries the Risch algorithm first, then falls back to rule-based methods if needed.

For using SymbolicIntegration.jl, see [the SymbolicIntegration.jl repository](https://github.com/JuliaSymbolics/SymbolicIntegration.jl) 
for more details. Quick examples:

```julia
using Symbolics
using SymbolicIntegration

@variables x

# Basic integrations (automatic method selection)
integrate(x^2, x)           # Returns (1//3)*(x^3)
integrate(1/x, x)           # Returns log(x)
integrate(exp(x), x)        # Returns exp(x)
integrate(1/(x^2 + 1), x)   # Returns atan(x)

# Rational function with complex roots (typically uses Risch)
f = (x^3 + x^2 + x + 2)/(x^4 + 3*x^2 + 2)
integrate(f, x)  # Returns (1//2)*log(2 + x^2) + atan(x)

# Nested transcendental functions (uses Risch)
integrate(1/(x*log(x)), x)  # Returns log(log(x))

# Explicit Risch method selection
integrate(sin(x)*cos(x), x, RischMethod())

# Configurable Risch method
risch = RischMethod(use_algebraic_closure=true, catch_errors=false)
integrate(exp(x)/(1 + exp(x)), x, risch)

# Rule-based method for cases not handled by Risch
integrate(x^(1/2), x)       # Uses rule-based method for non-integer powers
integrate(sqrt(x + 1), x)   # Handles algebraic functions via rules
```

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
