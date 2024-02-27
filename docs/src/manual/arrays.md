# [Symbolic Arrays](@id symbolic_arrays)

## Symbolic Arrays vs Arrays of Symbolic Expressions

Symbolics.jl contains two forms for handling symbolic arrays:

1. Arrays of symbolic expressions: these are Julia arrays with Symbolics.jl objects in them.
2. Symbolic Arrays: these are symbolic (O(1)) representations of arrays.

Arrays of symbolic expressions are simply Symbolics.jl objects put into Julia arrays. For
example:

```@example arrays
using Symbolics
@variables x y
u = [x,y]
```

is a vector of two symbolic variables. As shorthand,

```@example arrays
u2 = Symbolics.variables(:x, 1:3, 3:6)
```

creates a Julia matrix of symbolic variables. Indexing `u` or `u2` gives symbolic values
which act as a normal scalar symbolic value. This form these uses Julia's array functionality
and performs symbolic operations on the scalar values.

On the otherhand, Julia's symbolic array form is an O(1) representation of the whole array.

```@example arrays
@variables A[1:5, 1:3]
```

When using this form, `A[1,1]` is not a symbolic variable but a symbolic expression for
indexing the variable `A`. This representation holds linear algebra expressions in a
non-expanded form. For example:

```@example arrays
@variables B[1:3, 1:3]
A * B
```

in comparison to:

```@example arrays
a = Symbolics.variables(:a, 1:5, 1:3)
b = Symbolics.variables(:b, 1:3, 1:3)
a * b
```

This makes the symbolic array form much more efficient, but requires that the expressions
uses things with registered symbolic array functions which currently has much lower coverage.
Also, there are many fallbacks for which arrays of symbolics which makes this approach
more accessible but with larger expressions.

We recommend defaulting to arrays of symbolics unless you need the expression symplifications
of the symbolic array approach.

## Using Symbolic Arrays

Symbolic array-valued expressions (symbolic arrays) are supported by Symbolics. Symbolic array expressions propagate useful metadata that depends on input arrays: array dimension, element type and shape.

You can create a symbolic array variable with the following syntax:

```@example arrays
using Symbolics
@variables A[1:5, 1:3] b[1:3]
```

Here, `A` is a symbolic matrix of size `(5, 3)` and `b` is a symbolic vector of length 3.

```@example arrays
size(A)
```
```@example arrays
size(b)
```
```@example arrays
ndims(A)
```
```@example arrays
ndims(b)
```
```@example arrays
eltype(A)
```
```@example arrays
eltype(b)
```

## Array operations

Operations on symbolic arrays return symbolic array expressions:

```@example arrays
c = A * b
```
```@example arrays
size(c)
```
```@example arrays
eltype(c)
```

Adjoints, matrix-matrix, and matrix-vector multiplications are supported. Dot product returns a scalar-valued expression:

```@example arrays
b'b
```
```@example arrays
size(b'b)
```

Outer product returns a matrix:

```@example arrays
b * b'
```
```@example arrays
size(b*b')
```

### Broadcast, map and reduce


```@example arrays
A .* b'
```
```@example arrays
map(asin, (A*b))
```
```julia
#sum(A) #latexify not working
```
```@example arrays
typeof(sum(A))
```
```@example arrays
typeof(sum(A, dims=2))
```

### Indexing and delayed computation

Indexing array expressions is fairly flexible in Symbolics. Let's go through all the possible ways to index arrays.

#### Scalar indexing and scalarization

```@example arrays
AAt = A*A'
```
```@example arrays
AAt[2,3]
```

Here we indexed for the element (2,3), but we got back a symbolic indexing expression. You may want to force the element to be computed in terms of the elements of A. This can be done, using the `scalarize` function.

```@example arrays
Symbolics.scalarize(AAt[2,3])
```
```julia
@syms i::Int j::Int
Symbolics.scalarize(AAt[i,j])
```

In general, any scalar expression which is derived from array expressions can be scalarized.

```julia
#sum(A[:,1]) + sum(A[2,:])#latexify not working
```
```julia
Symbolics.scalarize(sum(A[:,1]) + sum(A[2,:]))
```
