# Expression Manipulation

Symbolics.jl provides functionality for easily manipulating expressions.
Most of the functionality comes by the expression objects obeying the standard
mathematical semantics. For example, if one has `A` a matrix of symbolic
expressions wrapped in `Num`, then `A^2` calculates the expressions for the
squared matrix.  It is thus encouraged to use standard Julia
for performing many of the manipulation on the IR. For example,
calculating the sparse form of the matrix via `sparse(A)` is valid, legible,
and easily understandable to all Julia programmers.

## Functionality Inherited From SymbolicUtils.jl

```@docs
SymbolicUtils.substitute
SymbolicUtils.simplify
```
Documentation for `rewriter` can be found [here](https://symbolicutils.juliasymbolics.org/rewrite/), using the `@rule` macro or the `@acrule` macro from SymbolicUtils.jl.

## Additional Manipulation Functions

Other additional manipulation functions are given below.

```@docs
Symbolics.get_variables
Symbolics.tosymbol
Symbolics.diff2term
Symbolics.degree
Symbolics.coeff
Symbolics.replace
Symbolics.occursin
Symbolics.filterchildren
Symbolics.fixpoint_sub
Symbolics.fast_substitute
Symbolics.symbolic_to_float
Symbolics.terms(x)
Symbolics.factors(x)
numerator(x::Union{Num, Symbolics.Symbolic})
denominator(x::Union{Num, Symbolics.Symbolic})
```
