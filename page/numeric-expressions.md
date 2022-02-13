@def title="Numeric expressions"

# Numeric expressions

The most basic symbolic expressions are symbolic variables. They can be created using the `@variables` macro.


\repl{
using Symbolics
@variables x y z
}

Arithmetic on variables returns symbolic expressions:

\repl{
x + y
sin(y) / sin(z)
}

**Simplifying assumptions**

By default, variables created with `@variables` macro make some _simplifying_ assumptions. Namely, multiplication and addition are assumed to be commutative and associative, divison assumes that the domain of the denominator does not contain zero.
These assumptions help Symbolics maintain a compact canonical representation of expressions. This saves on algorithmic complexity, memory and results in most expressions being simplified enough that calling `simplify` explicitly becomes unnecessary. (`simplify` is still necessary to apply things like trigonometric identities or to simplify fractions by eliminating the GCD from the numerator and denominator)

\repl{
x + y
y + x
x + y + y + x
x^2 * y * x
x/(x*y)
}

**Avoiding simplifying assumptions**

In some cases, you explicitly do not want these assumptions to be made. For example, if you wish to track the exact operations a symbolic expression was created with so as to apply some analysis that depends on the same, you may want to avoid the commutative and associative assumptions. In such cases, you can instantiate the variables of the `LiteralReal` type by specifying the same in `@variables` declaration:

\repl{
@variables a::LiteralReal b::LiteralReal c::LiteralReal
}

Arithmetic on these variables will not canonicalize the expressions:

\repl{
a + b
b + a
a + b + b + a
a^2 * b * a
a/a
a/(a*b)
}


## Accessing expression trees

The datastructure used to represent expressions may change to allow for faster canonicalization and efficient storage, and therefore, expressions need to be treated as black-box objects.

Another peculiarity is that every expression of the Real symbolic type is wrapped in a type called `Num` which is a subtype of `Base.Real`. This allows the use of these expressions in generic Julia functions meant for `Real` numbers. So to access the underlying expression tree, you must unwrap the expression first.

\repl{
ex = x + 2y
typeof(ex)
using Symbolics: unwrap
unwrap(ex)
typeof(unwrap(ex))
}

Notice that the expression tree of `x + 2y` is of the type `Add`. There are many such types to represent expressions. This is why Symbolics stipulates a standard set of _accessor functions_ to inspect expressions of every type uniformly.

Basically, every expression that is not a symbolic variable is an application of some Julia function on argument expressions. The accessor functions let us access the function and the arguments:

- `istree(ex)` tells you if an expression is a non-leaf node (i.e., not a variable or a number but something that is a function of further arguments).
- `operation(ex)` gives the julia function object that is being applied,
- `arguments(ex)` gives the arguments the function is being applied on.

\repl{
using Symbolics: istree, operation, arguments

istree(42) # a number is not a tree
istree(x)  # nor is a variable
ex = 2x + y
ex = unwrap(ex)
typeof(ex)
istree(ex)
operation(ex)
arguments(ex)
ex = 3x^2*y
ex = unwrap(ex)
typeof(ex)
istree(ex)
operation(ex)
arguments(ex)
}

## Complex expressions

Complex number expressions are represented using Julia's native `Complex` type where the real and imaginary parts are `Real` symbolic expressions we saw so far. (These are wrapped in `Num`, but the `Complex` object itself is not.)

\repl{
ex = a+im*b
typeof(ex)
real(ex)
imag(ex)
}

`unwrap` on a complex expression returns a special tree whose `operation` is `Complex` and the arguments are the real and imaginary parts. This makes generic tree-accessing code work with complex expressions.

\repl{
ex = unwrap(ex)
istree(ex)
operation(ex)
arguments(ex)
}
