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
and then we build the C source (the default is `expression=Val{true}`):

```@example converting_to_C
ccode = build_function(du, u, p, t, target=Symbolics.CTarget(), expression=Val{true})
```

`CTarget` returns C source as a `String`. It does not invoke a C compiler or
return a Julia callable: `expression=Val{false}` throws an error. Compile the
generated source yourself, load it with `Libdl`, and call it via `ccall`.

!!! note
    The following compile-and-call steps require a C compiler such as `gcc`.
    They are shown as non-executing documentation so the docs build does not
    depend on a compiler being present in the documentation environment.

```julia
using Libdl

# Write the generated C to a temporary shared library
libdir = mktempdir()
libpath = joinpath(libdir, "lotka." * dlext)
open(`gcc -fPIC -O3 -xc -shared -o $libpath -`, "w") do io
    print(io, ccode)
end

# Call the C function (default fname is :diffeqf)
du = rand(2); du2 = rand(2)
u = rand(2)
p = rand(4)
t = rand()
ccall((:diffeqf, libpath), Cvoid,
      (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64),
      du, u, p, t)
lotka_volterra!(du2, u, p, t)
du == du2 # true!
```
