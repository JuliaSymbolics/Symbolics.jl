# Symbolic arrays

Symbolic array-valued expressions (symbolic arrays) are supported by Symbolics. Symbolic array expressions propagate useful metadata that depends on input arrays: array dimension, element type and shape.

You can create a symbolic array variable with the following syntax:

```@example arrays
using Symbolics
@variables A[1:5, 1:3] b[1:3]
```

Here `A` is a symbolic matrix of size `(5, 3)` and `b` is a symbolic vector of length 3.

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
```@example arrays
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

Here we indexed for the element (2,3), but we got back a symbolic indexing expression. You may want to force the element to be computed in terms of the elements of A. This can be done, using `scalarize` function.

```@example arrays
Symbolics.scalarize(AAt[2,3])
```
```@example arrays
@syms i::Int j::Int
Symbolics.scalarize(AAt[i,j])
```

In general any scalar expression which is derived from array expressions can be scalarized.

```@example arrays
#sum(A[:,1]) + sum(A[2,:])#latexify not working
```
```@example arrays
Symbolics.scalarize(sum(A[:,1]) + sum(A[2,:]))
```
