# I/O, Saving, and Latex

Note that Julia's standard I/O functionality can be used to save
Symbolics expressions out to files. For example, here we will
generate an in-place version of `f` and save the anonymous function to
a `.jl` file:

```@example io
using Symbolics
@variables u[1:3]
function f(u)
  [u[1]-u[3],u[1]^2-u[2],u[3]+u[2]]
end
ex1, ex2 = build_function(f(u),u)
write("function.jl", string(ex2))
```

Now we can do something like:

```@example io
g = include("function.jl")
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

### Per-variable LaTeX names

Attach a `latexwrapper` to a variable to choose its LaTeX rendering
independently of its Julia name (similar to SymPy's `Symbol(r'\omega_{0}')`):

```julia
using Symbolics, Latexify

@variables x t
@variables w0 [latexwrapper = s -> raw"\omega_{0}"]
@variables vx(x, t) [latexwrapper = s -> "v_{x}"]

latexify(vx + w0^2)
# v_{x}\left( x, t \right) + \omega_{0}^{2}
```

The wrapper is a function from the symbol's string name to a LaTeX fragment.
That fragment is emitted verbatim, so underscores and carets are not escaped.
Without `latexwrapper`, multi-character names are wrapped in `\mathtt{...}`
by default.
