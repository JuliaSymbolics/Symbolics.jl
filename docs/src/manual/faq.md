# Frequently Asked Questions

## Limits of Symbolic Computation

### Transforming my function to a symbolic equation has failed. What do I do?

If you see the error:

```
ERROR: TypeError: non-boolean (Num) used in boolean context
```

this is likely coming from an algorithm which cannot be traced into a purely
symbolic algorithm. Many numerical solvers, for instance, have this property. It
shows up when you're doing something like if `x < tol`. If x is a number, then
this is true or false. If x is a symbol, then it's `x < tol`, so Julia just cannot
know how many iterations to do and throws an error.

This shows up in adaptive algorithms, for example:

```@example faq
function factorial(x)
  out = x
  while x > 1
    x -= 1
    out *= x
  end
  out
end
```

The number of iterations this algorithm runs for is dependent on the value of
`x`, and so there is no static representation of the algorithm. If `x` is 5,
then it's `out = x*(x-1)*(x-2)*(x-3)*(x-4)`, while if `x` is 3, then it's
`out = x*(x-1)*(x-2)`. It should thus be no surprise that:

```@example faq
using Symbolics
@variables x
try
    factorial(x)
catch e
    e
end
```

fails. It's not that there is anything wrong with this code, but it's not going
to work because fundamentally this is not a symbolically-representable algorithm.

The space of algorithms which can be turned into symbolic algorithms is what we
call quasi-static, that is, there is a way to represent the algorithm as static.
Loops are allowed, but the amount of loop iterations should not require that you
know the value of the symbol `x`. If the algorithm is quasi-static, then Symbolics.jl
tracing will produce the static form of the code, unrolling the operations, and
generating a flat representation of the algorithm.

#### What can be done?

If you need to represent this function `f` symbolically, then you'll need to make
sure it's not traced and instead is directly represented in the underlying
computational graph. Just like how `sqrt(x)` symbolically does not try to
represent the underlying algorithm, this must be done to your `f`. This is
done by doing `@register_symbolic f(x)`. If you have to define things like derivatives to
`f`, then [the function registration documentation](@ref function_registration).

## Equality and set membership tests
Comparing symbols with `==` produces a symbolic equality, not a `Bool`. To produce a `Bool`, call `isequal`.

To test if a symbol is part of a collection of symbols, i.e., a vector, either create a `Set` and use `in`, e.g.
```@example faq
try 
    x in [x]
catch e
    e
end
```
```@example faq
x in Set([x])
```
```@example faq
any(isequal(x), [x])
```

If `==` is used instead, you will receive `TypeError: non-boolean (Num) used in boolean context`. What this error 
is telling you is that the symbolic `x == y` expression is being used where a `Bool` is required, such as
`if x == y`, and since the symbolic expression is held lazily this will error because the appropriate branch cannot
be selected (since `x == y` is unknown for arbitrary symbolic values!). This is why the check `isequal(x,y)` is
required, since this is a non-lazy check of whether the symbol `x` is always equal to the symbol `y`, rather than
an expression of whether `x` and `y` currently have the same value.

## Understanding the Difference Between the Julia Variable and the Symbolic Variable

In the most basic usage of Symbolics, the name of the Julia variable
and the symbolic variable are the same. For example, when we do:

```@example faq
@variables a
```

the name of the symbolic variable is `a` and same with the Julia variable. However, we can
de-couple these by setting `a` to a new symbolic variable, for example:

```@example faq
b = only(@variables(a))
```

Now the Julia variable `b` refers to the variable named `a`. However, the downside of this current
approach is that it requires that the user writing the script knows the name `a` that they want to
place to the variable. But what if for example we needed to get the variable's name from a file?

To do this, one can interpolate a symbol into the `@variables` macro using `$`. For example:

```@example faq
a = :c
b = only(@variables($a))
```

In this example, `@variables($a)` created a variable named `c`, and set this variable to `b`.

## [Why does `A \ b` give `NaN` or make `simplify` extremely slow?](@id faq_symbolic_backslash)

For arrays of symbolic scalars (`Matrix{Num}`), `A \ b` uses an LU factorization
with soft pivoting (`sym_lu`): each pivot is chosen by fewest expression terms, not
by numeric magnitude. A pivot that is symbolically nonzero can still be numerically
zero after substitution, which produces nested `0/0`-style fractions and `NaN`.
Back-substitution also nests divisions, so the first unknown (`sol[1]`) is usually
the deepest expression and the slowest (or impossible) to `simplify`. Related
`DivideError`s during simplification of `1/0`-like forms are tracked separately
(see [issue 878](https://github.com/JuliaSymbolics/Symbolics.jl/issues/878)).

If you need a closed-form solution that stays valid for every nonsingular numeric
specialization of `A`, prefer the Laplace (cofactor) paths, which divide only by
`det(A)` (once, at the end) instead of nesting a division at every elimination step:

```julia
using Symbolics, LinearAlgebra
A = Symbolics.@variables(A[1:4, 1:4])[1] |> Symbolics.scalarize
b = Symbolics.@variables(b[1:4])[1] |> Symbolics.scalarize
sol = -(inv(A) * b)
# or Cramer:
# d = det(A)
# sol = [-det(hcat(A[:, 1:(i-1)], b, A[:, (i+1):end])) / d for i in 1:4]
```

Laplace expansion is factorial in the matrix size, so this route is practical only
for small dense systems. Timed on current Symbolics (Julia 1.12), a second call to
`inv(A)` for a fully symbolic `n×n` matrix is on the order of milliseconds through
`5×5`, about a second at `7×7`, about ten seconds at `8×8`, and about a minute and
a half at `9×9`; larger sizes grow quickly. Within the small-`n` range the expressions
evaluate correctly even when `A[1,1] == 0`, where `A \ b` can return `NaN`/`Inf`.