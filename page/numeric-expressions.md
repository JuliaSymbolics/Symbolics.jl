@def title="Numeric expressions"

# Expressions

The most basic symbolic expressions are symbolic variables. They can be created using the `@variables` macro.


\repl{
using Symbolics
@variables x y z;
x+y
}

Arithmetic on variables returns symbolic expressions:



Every expression in Symbolics has a symbolic type. While the expression itself is a rather opaque data structure, its symbolic type can be queried by calling the function [`symtype`](/api#symtype) on it. Here we refer to expressions whose symbolic type is a subtype of `Number` as **numeric expressions**.

The most basic type of expression is just a symbolic variable (sometimes just called a symbol). Symbolic variables can be created with the `@variables` macro. The macro takes a list of names and generates a vector of symbols with those names. It also creates Julia variables of the same names and assigns the symbolic variables to them correspondingly.

For example,
\repl{
1+1
2+3
}
```julia:repl1
using Symbolics

@variables w z α β
```

This creates variables named `w`, `z`, `α` and `β` and assigns them to a symbolic object with the same name. By default `@variables` creates variables of the symbolic type `Real`.

```julia:repl1
symtype(w)
```

Arithmetic on numeric expressions returns numeric expressions. To save space and time, Symbolics uses efficient datastructures to store expressions in a compact simplified form when they are created. This is illustrated in the output of the following simple arithmetic operations:

```julia:repl1
w + z + z

w + z * z

2(w + z)

w / w

w - w
```
\textoutput{repl1}

Note that the default regime of arithmetic on variables created by `@variables` makes several simplifying assumptions. They are,

# Safer Real expressions

Instead of `Real` variables, one can create `SafeReal` or `LiteralReal` variables which make lesser assumptions about their arithmatic. A `SafeReal` variable does not assume that division of expressions never results in a divison by zero.

## SafeReal
