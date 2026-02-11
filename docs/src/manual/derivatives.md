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
is_derivative
```

!!! note
    For symbolic differentiation, all registered functions in the symbolic expression
    need a registered derivative. For more information, see the
    [function registration](@ref function_registration) page.

## Substitution in Derivative Expressions

As of Symbolics.jl v7, [`substitute`](@ref) no longer recurses into the arguments of
`Differential` expressions. For example, `substitute(D(x), Dict(x => y))` will return
`D(x)` unchanged. This is because SymbolicUtils.jl v4 treats `Operator` subclasses
(including `Differential`) as substitution boundaries by default.

To substitute inside `Differential` applications, use [`substitute_in_deriv`](@ref) or
[`substitute_in_deriv_and_depvar`](@ref):

```julia
using Symbolics
@variables t x(t) y(t)
D = Differential(t)

substitute(D(x), Dict(x => y))                   # D(x) — does NOT penetrate
substitute_in_deriv(D(x), Dict(x => y))           # D(y) — penetrates Differential
substitute_in_deriv_and_depvar(D(x), Dict(x => y)) # D(y) — also penetrates dependent variables
```

```@docs
Symbolics.substitute_in_deriv
Symbolics.substitute_in_deriv_and_depvar
```

## High-Level Differentiation Functions

The following functions are not exported and thus must be accessed in a namespaced
way, i.e. `Symbolics.jacobian`.

```@docs
Symbolics.derivative
Symbolics.jacobian
Symbolics.sparsejacobian
Symbolics.sparsejacobian_vals
Symbolics.gradient
Symbolics.hessian
Symbolics.sparsehessian
Symbolics.sparsehessian_vals
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
function via the macro:

```julia
@register_derivative my_function(args...) i begin
    # ...
end
```

where `i` means that it's the derivative with respect to the `i`th argument. `args` is the
array of arguments, so, for example, if your function is `f(x,t)`, then `args = [x,t]`.
You should return an `Term` for the derivative of your function.

For example, `sin(t)`'s derivative (by `t`) is given by the following:

```@example derivatives
@register_derivative sin(t) 1 cos(t)
```

For more information see [`@register_derivative`](@ref).
