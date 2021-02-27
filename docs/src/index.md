# Symbolics.jl

Symbolics.jl is a fast and modern Computer Algebra System (CAS) for a fast and modern
programming language (Julia). The goal is to have a high-performance and parallelized
symbolic algebra system that is directly extendable in the same language as the users.

## Installation

To install Symbolics.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("Symbolics")
```

## Feature Summary

Because Symbolics.jl is built into the Julia language and works
with its dispatches, generic functions in Base Julia will work with symbolic
expressions! Make matrices of symbolic expressions and multiply them: it will
just work. Take the LU-factorization. Etc. Thus see
[the Julia Documentation](https://docs.julialang.org/en/v1/) for a large list
of functionality available in Symbolics.jl.

In addition, Symbolics.jl adds specific functions for:

- Simplification and substitution
- Differentiation
- Automatic conversion of Julia code to symbolic code
- Generation of (high performance and parallel) functions from symbolic expressions
- Automated sparsity detection of Julia code

and much more.

## Extension Packages

- [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl): Symbolic
  representations of models, like ODEs, SDEs, PDEs, optimization, optimal control,
  nonlinear solves, and more.
