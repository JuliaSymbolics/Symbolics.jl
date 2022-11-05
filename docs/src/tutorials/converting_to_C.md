# Automatic Conversion of Julia Code to C Functions

Since Symbolics.jl can trace Julia code into Symbolics IR that can be built and
compiled via `build_function` to C, this gives us a nifty way to automatically
generate C functions from Julia code! To see this in action, let's start with
[the Lotka-Volterra equations](https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations):

```@example converting_to_C
using Symbolics
function lotka_volterra!(du, u, p, t)
  x, y = u
  α, β, δ, γ = p
  du[1] = dx = α*x - β*x*y
  du[2] = dy = -δ*y + γ*x*y
end
```

Now we trace this into Symbolics:

```@example converting_to_C
@variables t du[1:2] u[1:2] p[1:4]
du = collect(du)
lotka_volterra!(du, u, p, t)
du
```
and then we build the function:

```@example converting_to_C
build_function(du, u, p, t, target=Symbolics.CTarget())
```

If we want to compile this, we do `expression=Val{false}`:

```@example converting_to_C
f = build_function(du, u, p, t, target=Symbolics.CTarget(), expression=Val{false})
```

now we check it computes the same thing:

```@example converting_to_C
du = rand(2); du2 = rand(2)
u = rand(2)
p = rand(4)
t = rand()
f(du, u, p, t)
lotka_volterra!(du2, u, p, t)
du == du2 # true!
```
