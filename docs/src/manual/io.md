# I/O, Saving, and Latex

Note that Julia's standard I/O functionality can be used to save
Symbolics expressions out to files. For example, here we will
generate an in-place version of `f` and save the anonymous function to
a `.jl` file:

```julia
using Symbolics
@variables u[1:3]
function f(u)
  [u[1]-u[3],u[1]^2-u[2],u[3]+u[2]]
end
ex1, ex2 = build_function(f(u),u)
write("function.jl", string(ex2))
```

Now we can do something like:

```julia
f = include("function.jl")
```

and that will load the function back in. Note that this can be done
to save the transformation results of Symbolics.jl so that
they can be stored and used in a precompiled Julia package.

## Latexification

Symbolics.jl's expressions support Latexify.jl, and thus

```julia
using Latexify
latexify(ex)
```

will produce LaTeX output from Symbolics models and expressions.
This works on basics like `Term` all the way to higher primitives
like `ODESystem` and `ReactionSystem`.
