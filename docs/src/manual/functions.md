# [Function Registration and Tracing](@id function_registration)

## Direct Tracing

Because Symbolics expressions respect Julia semantics, one way
to generate symbolic expressions is to simply place Symbolics
variables as inputs into existing Julia code. For example, the following
uses the standard Julia function for the Lorenz equations to generate
the symbolic expression for the Lorenz equations:

```@example registration
using Symbolics
function lorenz(du,u,p,t)
 du[1] = 10.0(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end
@variables t p[1:3] u(t)[1:3]
du = Array{Any}(undef, 3)
lorenz(du,u,p,t)
du
```
Or similarly:

```@example registration
@variables t x(t) y(t) z(t) dx(t) dy(t) dz(t) σ ρ β
du = [dx,dy,dz]
u = [x,y,z]
p = [σ,ρ,β]
lorenz(du,u,p,t)
du
```

## Registering Functions

The Symbolics graph only allows registered Julia functions within its type.
All other functions are automatically traced down to registered
functions. By default, Symbolics.jl pre-registers the common functions
utilized in [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl)
and pre-defines their derivatives. However, the user can utilize the
[`@register_symbolic`](@ref) macro to add their function to allowed functions
of the computation graph.

Additionally, [`@register_array_symbolic`](@ref) can be used to define array functions.
For size propagation it's required that a computation of how the sizes are computed is
also supplied.

## Defining Derivatives of Registered Functions

In order for symbolic differentiation to work, an overload of `Symbolics.derivative` is
required. The syntax is `derivative(typeof(f), args::NTuple{i,Any}, ::Val{j})` where
`i` is the number of arguments to the function and `j` is which argument is being
differentiated. So for example:

```julia
function derivative(::typeof(min), args::NTuple{2,Any}, ::Val{1})
    x, y = args
    IfElse.ifelse(x < y, one(x), zero(x))
end
```

is the partial derivative of the Julia `min(x,y)` function with respect to `x`.

!!! note
    Downstream symbolic derivative functionality only work if every partial derivative that
    is required in the derivative expression is defined. Note that you only need to define
    the partial derivatives which are actually computed.

## Registration API

```@docs
@register_symbolic
@register_array_symbolic
```
