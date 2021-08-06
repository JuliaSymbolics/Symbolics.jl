# Automatic Conversion of Julia Code to C Functions

Since Symbolics.jl can trace Julia code into Symbolics IR that can be built and
compiled via `build_function` to C, this gives us a nifty way to automatically
generate C functions from Julia code! To see this in action, let's start with
[the Lotka-Volterra equations](https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations):

```julia
using Symbolics
function lotka_volterra!(du, u, p, t)
  x, y = u
  α, β, δ, γ = p
  du[1] = dx = α*x - β*x*y
  du[2] = dy = -δ*y + γ*x*y
end
```

Now we trace this into Symbolics:

```julia
@variables t du[1:2] u[1:2] p[1:4]
du = collect(du)
lotka_volterra!(du, u, p, t)
```

which gives:

```julia
du = Num[p₁ * u₁ - (p₂ * u₁) * u₂, -p₃ * u₂ + (p₄ * u₁) * u₂]
```

and then we build the function:

```julia
build_function(du, u, p, t, target=Symbolics.CTarget())

void diffeqf(double* du, double* RHS1, double* RHS2, double RHS3) {
  du[0] = RHS2[0] * RHS1[0] - (RHS2[1] * RHS1[0]) * RHS1[1];
  du[1] = -(RHS2[2]) * RHS1[1] + (RHS2[3] * RHS1[0]) * RHS1[1];
}
```

If we want to compile this, we do `expression=Val{false}`:

```julia
f = build_function(du, u, p, t, target=Symbolics.CTarget(), expression=Val{false})
```

now we check it computes the same thing:

```julia
du = rand(2); du2 = rand(2)
u = rand(2)
p = rand(4)
t = rand()
f(du, u, p, t)
lotka_volterra!(du2, u, p, t)
du == du2 # true!
```
