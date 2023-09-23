# Supported types and dispatch in Symbolics

There is a tension between types as a representation of expression trees and types that are a subtype of types already present in Julia.

We want to be able to deal with expression trees in a unified way and not constrain expression trees themselves to be under an abstract type in Julia's type hierarchy. (For example, if we said that all expression trees are subtype of Real, then we couldn't represent array operations using the same expression tree.). But we also want to be able to pass in expression trees into places in existing code that accept `Real` values.

We accomplish this by **wrapping** expression trees in a simple wrapper type which is a subtype of our desired abstract type. For example, we wrap expression trees in the type `Num` which is a subtype of `Real` to make it behave like a Real number.

The methods on `Num` objects are forwarded to the wrapped expression tree. And care is taken so that an expression tree never internally contains `Num` -- this is both for performance and separation of concerns.

User-facing APIs in Symbolics always take wrapped objects like `Num`, they are then internally unwrapped for expression tree manipulation.

Due to it requiring such wrappers, we only fully support a limited number of types as both the types of expression trees and the type as Julia sees them.

These types are

- Real numbers (wrapped using `Num`)
- complex numbers (stored as `Complex{Num}` where `Complex` is from Base Julia)
- arrays of Real and complex numbers (wrapped using `Arr`, so `Arr{Num}` or `Arr{Complex{Num}}`)

## `@variables` and types

Use the syntax `@variables x::T` to create a symbol named `x` of symbolic type `T`. If `T` is a subtype of any of the above listed types which support a wrapper, the resulting variable will be wrapped in that type. As seen in the examples below, x,z,X,Z all have a suitable wrapper type. Hence, their types are shown. However, `s` being of symbolic type `String` does not have a corresponding wrapper supported by Symbolics, and hence, it returns a `Sym{String}` object. This is the trivial expression tree of a single variable without a wrapper, and is not a subtype of String or AbstractString.

```@example types
using Symbolics
@variables x::Real z::Complex{Real} (X::Real)[1:10, 1:10] (Z::Complex{Real})[1:10] s::String
```
```@example types
typeof(x)
```
```@example types
typeof(z)
```
```@example types
typeof(X)
```
```@example types
typeof(Z)
```
```@example types
typeof(s)
```

## Canonicalization, and avoiding it

By default, Symbolics will simplify expressions to an internal canonical form upon construction:

```@example types
@variables x y
(x + y, y + x) # assumes commutativity of +
```
```@examples types
x + y + y + x # assumes associativity also
```
```@examples types
x^2 * y * x # same goes for *
```
```@examples types
x/(x*y) # assumes x!=0
```

These assumptions are made so that the obvious simplifications are done as soon as someone writes them. If you want to ask Symbolics to not make these assumptions, you can annotate them with the `::LiteralReal` type when declaring them with `@variables`. By default, variables are implicitly annotated by `::Real`, `::LiteralReal` is an internal type which is a subtype of `Real` which is used to create expressions with the property that they are not auto-simplified.

```@examples types
@variables a::LiteralReal b::LiteralReal

(a + b,
 b + a,
 a + b + b + a,
 a^2 * b * a,
 a/(a*b))

You can later force expressions to be simplified by using the `simplify` function.


```@examples types
simplify(a + b + b + a)
```
