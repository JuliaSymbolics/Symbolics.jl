# Symbolic Limits

Experimental symbolic limit support is provided by the [`limit`](@ref) function, documented
below. See [SymbolicLimits.jl](https://github.com/SciML/SymbolicLimits.jl) for more
information and implementation details.

```@docs
limit
```

### SymPy Integration

SymPy also includes limits as well, and the SymPy.jl extensions allow for automatically converting
Symbolics expressions for use in its solvers.

```@docs
Symbolics.sympy_limit
```