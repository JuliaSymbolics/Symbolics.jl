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
NormalDerivative
expand_derivatives
is_derivative
```

!!! note

    For symbolic differentiation, all registered functions in the symbolic expression
    need a registered derivative. For more information, see the
    [function registration](@ref function_registration) page.

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
Symbolics.get_differential_vars
Symbolics.substitute_in_deriv
Symbolics.substitute_in_deriv_and_depvar
```

## Vector Calculus Operators

Symbolics also provides multivariate differential operators for vector
calculus: `Grad` (`∇`), `Divergence` (`∇⋅`),
`Curl` (`∇×`) and `Laplacian` (`∇⋅∇`). Each is parameterized by
`E`: `Op{false}` (the default, also constructed by `Op(vars)`) builds a lazy
term with the operator as the head, which `expand_derivatives` expands, while
`Op{true}` expands into `Differential`s eagerly at application. `Grad`
composes with `LinearAlgebra.dot` and `LinearAlgebra.cross` into the other
operators. `Curl` accepts two or three independent variables: on a
2-component field it produces the scalar curl `D(vars[1])(v[2]) - D(vars[2])(v[1])`, and on a 3-component field the usual vector curl:

```julia
julia> using Symbolics, LinearAlgebra

julia> @variables x y z u[1:3];

julia> ∇ = Grad([x, y, z]);

julia> ∇(x^2 + y^2 + z^2)          # gradient, lazily
(∇(x^2 + y^2 + z^2))[1:3]

julia> expand_derivatives(ans)     # expanded
(SymbolicUtils.array_literal((3,), 2x, 2y, 2z))[1:3]

julia> ∇ ⋅ u                        # divergence
∇⋅(u)

julia> ∇ × u                        # curl
(∇×(u))[1:3]

julia> (∇ ⋅ ∇)(x^2 + y^2 + z^2)     # Laplacian
∇⋅(∇(x^2 + y^2 + z^2))
```

With two independent variables, `∇ × v` is the scalar curl of a 2-vector field,
identical to `Curl([x, y])(v)`:

```julia
julia> expand_derivatives(Grad([x, y]) × [x^2 * y, y^3 + x])
1 - (x^2)
```

```@docs
Symbolics.Grad
Symbolics.Divergence
Symbolics.Curl
Symbolics.Laplacian
```

## Adding Analytical Derivatives

There are many derivatives pre-defined by
[DiffRules.jl](https://github.com/JuliaDiff/DiffRules.jl).
For example,

```@example derivatives
using Symbolics
@variables x y z
f(x, y, z) = x^2 + sin(x+y) - z
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

## Declaring Total Derivatives

An opaque function can declare a total derivative independently of its partial
derivatives by specializing [`Symbolics.total_derivative`](@ref) on its own function
type. For example, a position function can have partial derivatives `1` and `0`,
while declaring velocity as its time derivative:

```@example derivatives
@variables time s(time) v(time)
pos(s, v) = s
@register_symbolic pos(s, v)
@register_derivative pos(s, v) 1 Symbolics.SConst(1)
@register_derivative pos(s, v) 2 Symbolics.SConst(0)

function Symbolics.total_derivative(::typeof(pos), args, D)
    isequal(D.x, Symbolics.unwrap(time)) || return nothing
    return args[2]
end

Dt = Differential(time)
expand_derivatives(Dt(pos(s, v))) # v
Symbolics.jacobian([pos(s, v)], [s, v]) # [1 0]
```

The hook receives read-only, unwrapped arguments and a first-order `Differential`.
Return `nothing` for other differentiation variables to retain the usual chain
rule and registered partial derivatives. The hook also works inside products and
nonlinear expressions. Higher-order and mixed derivatives differentiate the returned
expression normally; the declared derivative need not equal the chain rule applied
to the registered partials.

Every `Differential(time)` uses the hook, including uses intended as a partial
derivative with respect to `time`. For example,
`Symbolics.jacobian([pos(x, y)], [time])` returns the declared value `[y;;]`,
and ModelingToolkit's `tgrad` would use it too. Define a hook only where this is
the intended semantics. Mixed derivatives need not commute when a hook is used.

A hook must not return `D(f(args...))` for the same call: expanding that result
invokes the hook again and recurses forever.

```@docs
Symbolics.total_derivative
```
