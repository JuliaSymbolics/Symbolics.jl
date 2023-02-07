# Derivatives and Differentials

A `Differential(op)` is a partial derivative with respect to `op`,
which can then be applied to some other operations. For example, `D=Differential(t)`
is what would commonly be referred to as `d/dt`, which can then be applied to
other operations using its function call, so `D(x+y)` is `d(x+y)/dt`.

By default, the derivatives are left unexpanded to capture the symbolic
representation of the differential equation. If the user would like to expand
out all the differentials, the `expand_derivatives` function eliminates all
the differentials down to basic one-variable expressions.

```@docs
Differential
expand_derivatives
```

## High-Level Differentiation Functions

The following functions are not exported and thus must be accessed in a namespaced
way, i.e. `Symbolics.jacobian`.

```@docs
Symbolics.derivative
Symbolics.jacobian
Symbolics.sparsejacobian
Symbolics.gradient
Symbolics.hessian
Symbolics.sparsehessian
```

## Adding Analytical Derivatives

There are many derivatives pre-defined by
[DiffRules.jl](https://github.com/JuliaDiff/DiffRules.jl).
For example,
```@example derivatives
using Symbolics
@variables x y z
f(x,y,z) = x^2 + sin(x+y) - z
```

`f` automatically has the derivatives defined via the tracing mechanism. It will do
this by directly building the internals of your function and
differentiating that.

However, often you may want to define your own derivatives so that way
automatic Jacobian etc. calculations can utilize this information. This can
allow for more succinct versions of the derivatives to be calculated
to scale to larger systems. You can define derivatives for your
function via the dispatch:

```julia
# `N` arguments are accepted by the relevant method of `my_function`
Symbolics.derivative(::typeof(my_function), args::NTuple{N,Any}, ::Val{i})
```

where `i` means that it's the derivative with respect to the `i`th argument. `args` is the
array of arguments, so, for example, if your function is `f(x,t)`, then `args = [x,t]`.
You should return an `Term` for the derivative of your function.

For example, `sin(t)`'s derivative (by `t`) is given by the following:

```@example derivatives
Symbolics.derivative(::typeof(sin), args::NTuple{1,Any}, ::Val{1}) = cos(args[1])
```
