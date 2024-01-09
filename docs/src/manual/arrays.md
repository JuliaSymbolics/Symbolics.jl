# Symbolic arrays

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

Tensor or Einstein notation can also be used for array operations using the `@arrayop` macro. For example, a matrix-vector product can be written as:

```julia
julia> d = Symbolics.@arrayop (i,) A[i,j]*b[j]
@arrayop(_[i] := A[i, j]*b[j])
```
These array operations are computed lazily:

```julia
julia> c = A*b
(A*b)[1:5]

julia> isequal(collect(c), collect(d))
true
```

This can be useful for vector calculus operations such as the gradient of a vector expression:

```julia
julia> @variables x y;

julia> V = [sin(x), x^2 + y]
2-element Vector{Num}:
  sin(x)
 y + x^2

julia> D = Differential.([x,y]);
julia> expand_derivatives.(collect(Symbolics.@arrayop (i,j) D[i](V[j])))
2×2 Matrix{Num}:
 cos(x)  2x
      0   1
```
or the divergence of a vector expression:
```julia
julia> expand_derivatives.(collect(Symbolics.@arrayop () D[i](V[i])))
1 + cos(x)
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

Here we indexed for the element (2,3), but we got back a symbolic indexing expression. You may want to force the element to be computed in terms of the elements of A. This can be done, using the `scalarize` function.

```@example arrays
Symbolics.scalarize(AAt[2,3])
```
```@example arrays
@syms i::Int j::Int
Symbolics.scalarize(AAt[i,j])
```

In general, any scalar expression which is derived from array expressions can be scalarized.

```@example arrays
#sum(A[:,1]) + sum(A[2,:])#latexify not working
```
```@example arrays
Symbolics.scalarize(sum(A[:,1]) + sum(A[2,:]))
```
