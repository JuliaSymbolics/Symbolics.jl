# Working with External Symbolics Packages: SymPy, Mathematica, Oscar, and Beyond

Symbolics.jl takes an inclusive development philosophy: yes there are other computer algebra systems
such as SymPy and Mathematica, yes they can do some things that Symbolics.jl cannot natively do, and
no that is not a bad thing. Instead, we maintain bridges to these other languages in order to encourage
users to take full advantage of other tools as part of their workflow!

The following are the wrappers and interactions users should be aware of.

## Using Symbolics.jl with SymPy 

SymPy is a Python-based symbolic library. While it is generally not computationally performant, it is
good in its completeness of features and as such it can be helpful to call its functions from time to
time. Because of this, Symbolics.jl maintains extensions for SymPy.jl which facilitates translations.
The extension is triggered by adding the SymPy.jl packages doing `using SymPy`. When this is done,
Symbolics.jl gives functions which allow for translating to and from SymPy expressions:

```@docs
Symbolics.symbolics_to_sympy
Symbolics.sympy_to_symbolics
```

In addition, many of the features pages include docstrings for functionality given by this bridge.
These functions all have `sympy` in the name, and includes:

* `sympy_simplify`
* `sympy_linear_solve`
* `sympy_algebraic_solve`
* `sympy_integrate`
* `sympy_limit`
* `sympy_simplify`
* `sympy_ode_solve`

## Using Symbolics.jl with Mathematica / Wolfram's MathLink

The library [SymbolicsMathLink.jl](https://github.com/eswagel/SymbolicsMathLink.jl) is capable of
converting Symbolics.jl expressions for use with Mathematica calls via Wolfram's MathLink.
The package requires an installation of either [Mathematica](http://www.wolfram.com/mathematica/) or the free [Wolfram Engine](https://www.wolfram.com/engine/). It will attempt to find the installation at build time; if this fails, please see the [installation troubleshoot on the MathLink.jl README](https://github.com/JuliaInterop/MathLink.jl?tab=readme-ov-file#installation-troubleshoot).

Round-trip conversion is given by the functions:

```julia
expr_to_mathematica(juliaSymbolicExpression)
mathematica_to_expr(W`Some Mathematica expression`)
```

Some examples in action include using Mathematica for equation solving:

```julia
julia> using SymbolicsMathLink

julia> @variables x;
julia> expr = x^2 + x - 1;

julia> result = wcall("Solve", expr~0)
2-element Array{Num,1}:
    -1 + x
    1 + x
```

and solving differential equations:

```julia
julia> @variables vars(x)[1:2];
julia> expr = Differential(x)(vars[1]) + 2
2 + Differential(x)((vars(x))[1])
julia> result = wcall("DSolveValue", expr~0, vars[1], x)
C_1 - 2x
```

## Groebner Basis via Groebner.jl

[Groebner.jl](https://github.com/sumiya11/Groebner.jl) is a very efficient system for computing
Groebner basis written in Julia. Its publication 
[details how it is much more efficient than Maple, msolve, Mathematica, and Singular](https://arxiv.org/pdf/2304.06935) for these computing a Groebner basis. As such, it is highly recommended that one use this
system for the computation. Symbolics.jl provides an extension which automatically converts Symbolics.jl
polynomials into its specialized form for the computation and converts the result back, and thus no work
is required on the users end for integrating.

```@docs
Symbolics.groebner_basis
```

## Nemo.jl

[Nemo.jl's abstract algebra functionality (and by extension FLINT)](https://nemocas.github.io/Nemo.jl/stable/) 
is used as part of many of the symbolic solver routines. If Nemo is required an error message specifying
the requirement to load Nemo.jl will be given. For more information, see the 
[symbolic solver page](@ref symbolic_solver)