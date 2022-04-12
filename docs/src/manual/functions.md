# [Function Registration and Tracing](@id function_registration)

## Direct Tracing

Because Symbolics expressions respect Julia semantics, one way
to generate symbolic expressions is to simply place Symbolics
variables as inputs into existing Julia code. For example, the following
uses the standard Julia function for the Lorenz equations to generate
the symbolic expression for the Lorenz equations:

```julia
function lorenz(du,u,p,t)
 du[1] = 10.0(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end
@variables t p[1:3] u[1:3](t) du[1:3](t)
lorenz(du,u,p,t)
du
```

```julia
3-element Array{Num,1}:
                 10.0 * (u₂(t) - u₁(t))
         u₁(t) * (28.0 - u₃(t)) - u₂(t)
u₁(t) * u₂(t) - 2.6666666666666665 * u₃(t)
```

Or similarly:

```julia
@variables t x(t) y(t) z(t) dx(t) dy(t) dz(t) σ ρ β
du = [dx,dy,dz]
u = [x,y,z]
p = [σ,ρ,β]
lorenz(du,u,p,t)
du
```

```julia
3-element Array{Num,1}:
                10.0 * (y(t) - x(t))
         x(t) * (28.0 - z(t)) - y(t)
x(t) * y(t) - 2.6666666666666665 * z(t)
```

## Registering Functions

The Symbolics graph only allows registered Julia functions within its type.
All other functions are automatically traced down to registered
functions. By default, Symbolics.jl pre-registers the common functions
utilized in [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl)
and pre-defines their derivatives. However, the user can utilize the
[`@register_symbolic`](@ref) macro to add their function to allowed functions
of the computation graph.

```@docs
@register_symbolic
```
