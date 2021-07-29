# Symbolics.jl

[![Github Action CI](https://github.com/JuliaSymbolics/Symbolics.jl/workflows/CI/badge.svg)](https://github.com/JuliaSymbolics/Symbolics.jl/actions)
[![Coverage Status](https://coveralls.io/repos/github/JuliaSymbolics/ModelingToolkit.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaSymbolics/Symbolics.jl?branch=master)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://symbolics.juliasymbolics.org/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://symbolics.juliasymbolics.org/dev/)

Symbolics.jl is a fast and modern Computer Algebra System (CAS) for a fast and modern
programming language (Julia). The goal is to have a high-performance and parallelized
symbolic algebra system that is directly extendable in the same language as the users.

## Installation

To install Symbolics.jl, use the Julia package manager:

```julia
julia> using Pkg
julia> Pkg.add("Symbolics")
```

## Documentation

For information on using the package, see the [stable documentation](https://juliasymbolics.github.io/Symbolics.jl/dev/).
Use the [in-development documentation](https://juliasymbolics.github.io/Symbolics.jl/dev/)
for the version of the documentation which contains the unreleased features.

## Relationship to Other Packages

- [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl): This is a
  rule-rewriting system that is the core of Symbolics.jl. Symbolics.jl builds off of
  SymbolicUtils.jl to extend it to a whole symbolic algebra system, complete with
  support for differentation, solving symbolic systems of equations, etc. If you're
  looking for the barebones to build a new CAS for specific algebras, SymbolicUtils.jl
  is that foundation. Otherwise, Symbolics.jl is for you.
- [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl): This is a
  symbolic-numeric modeling system for the SciML ecosystem. It heavily uses Symbolics.jl
  for its representation of symbolic equations along with tools like differentiation,
  and adds the representation of common modeling systems like ODEs, SDEs, and more.

## Example

```julia
julia> using Symbolics

julia> @variables t x y
julia> D = Differential(t)

julia> z = t + t^2
julia> D(z) # symbolic representation of derivative(t + t^2, t)
Differential(t)(t + t^2)

julia> expand_derivatives(D(z))
1 + 2t

julia> Symbolics.jacobian([x + x*y, x^2 + y],[x, y])
2×2 Matrix{Num}:
 1 + y  x
    2x  1

julia> B = simplify.([t^2 + t + t^2  2t + 4t
                  x + y + y + 2t  x^2 - x^2 + y^2])
2×2 Matrix{Num}:
  t + 2(t^2)   6t
 x + 2t + 2y  y^2

julia> simplify.(substitute.(B, (Dict(x => y^2),)))
2×2 Matrix{Num}:
    t + 2(t^2)   6t
 2t + y^2 + 2y  y^2

julia> substitute.(B, (Dict(x => 2.0, y => 3.0, t => 4.0),))
2×2 Matrix{Num}:
 36.0  24.0
 16.0   9.0
```

## Citation

If you use Symbolics.jl, please [cite this paper](https://arxiv.org/abs/2105.03949)

```bib
@article{gowda2021high,
  title={High-performance symbolic-numerics via multiple dispatch},
  author={Gowda, Shashi and Ma, Yingbo and Cheli, Alessandro and Gwozdz, Maja and Shah, Viral B and Edelman, Alan and Rackauckas, Christopher},
  journal={arXiv preprint arXiv:2105.03949},
  year={2021}
}
```
