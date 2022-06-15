# Symbolics.jl

Symbolics.jl is a fast and modern Computer Algebra System (CAS) for a fast and modern
programming language (Julia). The goal is to have a high-performance and parallelized
symbolic algebra system that is directly extendable in the same language as the user's.

## Installation

To install Symbolics.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("Symbolics")
```

## Citation

If you use Symbolics.jl, please [cite this paper](https://arxiv.org/abs/2105.03949)

```
@article{gowda2021high,
  title={High-performance symbolic-numerics via multiple dispatch},
  author={Gowda, Shashi and Ma, Yingbo and Cheli, Alessandro and Gwozdz, Maja and Shah, Viral B and Edelman, Alan and Rackauckas, Christopher},
  journal={arXiv preprint arXiv:2105.03949},
  year={2021}
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

## Extension Packages

Below is a list of known extension packages. If you would like your package
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
