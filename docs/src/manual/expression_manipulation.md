# Expression Manipulation

Symbolics.jl provides functionality for easily manipulating expressions.
Most of the functionality comes by the expression objects obeying the standard
mathematical semantics. For example, if one has `A` a matrix of symbolic
expressions wrapped in `Num`, then `A^2` calculates the expressions for the
squared matrix.  In that sense, it is encouraged that one uses standard Julia
for performing a lot of the manipulation on the IR, as, for example,
calculating the sparse form of the matrix via `sparse(A)` is valid, legible,
and easily understandable to all Julia programmers.

Other additional manipulation functions are given below.

```@docs
get_variables
substitute
tovar
toparam
tosymbol
makesym
diff2term
solve_for
```
