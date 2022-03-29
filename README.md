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
  support for differentiation, solving symbolic systems of equations, etc. If you're
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

If you use Symbolics.jl, please [cite this paper](https://dl.acm.org/doi/10.1145/3511528.3511535) 
(or see the free [arxiv version](https://arxiv.org/abs/2105.03949))

```bib
@article{10.1145/3511528.3511535,
author = {Gowda, Shashi and Ma, Yingbo and Cheli, Alessandro and Gw\'{o}\'{z}zd\'{z}, Maja and Shah, Viral B. and Edelman, Alan and Rackauckas, Christopher},
title = {High-Performance Symbolic-Numerics via Multiple Dispatch},
year = {2022},
issue_date = {September 2021},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
volume = {55},
number = {3},
issn = {1932-2240},
url = {https://doi.org/10.1145/3511528.3511535},
doi = {10.1145/3511528.3511535},
abstract = {As mathematical computing becomes more democratized in high-level languages, high-performance symbolic-numeric systems are necessary for domain scientists and engineers to get the best performance out of their machine without deep knowledge of code optimization. Naturally, users need different term types either to have different algebraic properties for them, or to use efficient data structures. To this end, we developed Symbolics.jl, an extendable symbolic system which uses dynamic multiple dispatch to change behavior depending on the domain needs. In this work we detail an underlying abstract term interface which allows for speed without sacrificing generality. We show that by formalizing a generic API on actions independent of implementation, we can retroactively add optimized data structures to our system without changing the pre-existing term rewriters. We showcase how this can be used to optimize term construction and give a 113x acceleration on general symbolic transformations. Further, we show that such a generic API allows for complementary term-rewriting implementations. Exploiting this feature, we demonstrate the ability to swap between classical term-rewriting simplifiers and e-graph-based term-rewriting simplifiers. We illustrate how this symbolic system improves numerical computing tasks by showcasing an e-graph ruleset which minimizes the number of CPU cycles during expression evaluation, and demonstrate how it simplifies a real-world reaction-network simulation to halve the runtime. Additionally, we show a reaction-diffusion partial differential equation solver which is able to be automatically converted into symbolic expressions via multiple dispatch tracing, which is subsequently accelerated and parallelized to give a 157x simulation speedup. Together, this presents Symbolics.jl as a next-generation symbolic-numeric computing environment geared towards modeling and simulation.},
journal = {ACM Commun. Comput. Algebra},
month = {jan},
pages = {92–96},
numpages = {5}
}
```
