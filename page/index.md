@def title = "Symbolics.jl — Fast, modern, symbolic-numeric programming"
@def hasmath = false
@def hascode = true
<!-- Note: by default hasmath == true and hascode == false. You can change this in
the config file by setting hasmath = false for instance and just setting it to true
where appropriate -->

~~~
<h1>Symbolics.jl</h1>
~~~

~~~
<p style="font-size: 1.25em; line-height: 1.67em; text-align: center; margin: 1em 0; color: #111;">
<a href="https://github.com/JuliaSymbolics/Symbolics.jl">Symbolics.jl</a> is a fast and modern Computer Algebra System (CAS) for a fast and modern programming language.
</p>
~~~

## Installation

To install Symbolics.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("Symbolics")
```

## Citation

If you use Symbolics.jl, please [cite this paper](https://dl.acm.org/doi/10.1145/3511528.3511535) 
(or see the free [arxiv version](https://arxiv.org/abs/2105.03949))

```

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

## Feature Summary

Because Symbolics.jl is built into the Julia language and works
with its dispatches, generic functions in Base Julia will work with symbolic
expressions! Make matrices of symbolic expressions and multiply them: it will
just work. Take the LU-factorization. Etc. Thus see
[the Julia Documentation](https://docs.julialang.org/en/v1/) for a large list
of functionality available in Symbolics.jl.

A general list of the features is:

- Symbolic arithmetic with type information and multiple dispatch
- Symbolic polynomials and trigonometric functions
- Pattern matching, simplification and substitution
- Differentiation
- Symbolic linear algebra (factorizations, inversion, determinants, eigencomputations, etc.)
- Discrete math (representations of summations, products, binomial coefficients, etc.)
- Logical and Boolean expressions
- Symbolic equation solving and conversion to arbitrary precision
- Support for non-standard algebras (non-commutative symbols and customizable rulesets)
- Special functions (list provided by [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl))
- Automatic conversion of Julia code to symbolic code
- Generation of (high performance and parallel) functions from symbolic expressions
- Fast automated sparsity detection and generation of sparse Jacobian and Hessians

and much more.

## Packages that build on Symbolics

Below is a list of known packages that build on Symbolics. If you would like your package
to be listed here, feel free to open a pull request!

- [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl): Symbolic representations of common numerical systems
    - Ordinary differential equations
    - Stochastic differential equations
    - Partial differential equations
    - Nonlinear systems
    - Optimization problems
    - Optimal Control
    - Causal and acausal modeling (Simulink/Modelica)
    - Automated model transformation, simplification, and composition
- [Catalyst.jl](https://github.com/SciML/Catalyst.jl): Symbolic representations of chemical reactions
    - Symbolically build and represent large systems of chemical reactions
    - Generate code for ODEs, SDEs, continuous-time Markov Chains, and more
    - Simulate the models using the SciML ecosystem with O(1) Gillespie methods
- [DataDrivenDiffEq.jl](https://github.com/SciML/DataDrivenDiffEq.jl): Automatic identification of equations from data
    - Automated construction of ODEs and DAEs from data
    - Representations of Koopman operators and Dynamic Mode Decomposition (DMD)
- [SymbolicRegression.jl](https://github.com/MilesCranmer/SymbolicRegression.jl): Distributed High-Performance symbolic regression
    - Parallelized generic algorithms for finding equations from data
    - Pareto frontier based scoring
- [ReversePropagation.jl](https://github.com/dpsanders/ReversePropagation.jl): Source-to-source reverse mode automatic differentiation
    - Automated tracing of code and construction of backpropagation equations
    - Composes with symbolic transformation and simplification functionality
