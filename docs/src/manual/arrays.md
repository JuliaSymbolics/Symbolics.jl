# Symbolic arrays

Symbolic array-valued expressions (symbolic arrays) are supported by Symbolics. Symbolic array expressions propagate useful metadata that depends on input arrays: array dimension, element type and shape.

You can create a symbolic array variable with the following syntax:

```julia
@variables A[1:5, 1:3] b[1:3]
```

Here `A` is a symbolic matrix of size `(5, 3)` and `b` is a symbolic vector of length 3.

```julia
julia> size(A)
(5, 3)

julia> size(b)
(3,)

julia> ndims(A)
2

julia> ndims(b)
1

julia> eltype(A)
Real

julia> eltype(b)
Real
```

## Array operations

Operations on symbolic arrays return symbolic array expressions:

```julia
julia> c = A * b
(A*b)[1:5]

julia> size(c)
(5,)

julia> eltype(c)
Real
```

Adjoints, matrix-matrix, and matrix-vector multiplications are supported. Dot product returns a scalar-valued expression:

```julia
julia> b'b
adjoint(b)*b

julia> size(b'b)
()
```

Outer product returns a matrix:

```julia
julia> b * b'
(b*adjoint(b))[1:3,1:3]

julia> size(b*b')
(3, 3)
```

### Broadcast, map and reduce


```julia
julia> A .* b'
(broadcast(*, A, adjoint(b)))[1:5,1:3]
```

```julia
julia> map(asin, (A*b))
(map(asin, A*b))[1:5]
```

```julia
julia> sum(A)
...

julia> typeof(sum(A))
Num # it's a scalar!

julia> typeof(sum(A, dims=2))
Arr{Real, 2} # it's a vector
```

### Indexing and delayed computation

Indexing array expressions is fairly flexible in Symbolics. Let's go through all the possible ways to index arrays.

#### Scalar indexing and scalarization

```julia
julia> AAt = A*A'
(A*adjoint(A))[1:5,1:5]

julia> AAt[2,3]
(A*adjoint(A))[2,3]
```

Here we indexed for the element (2,3), but we got back a symbolic indexing expression. You may want to force the element to be computed in terms of the elements of A. This can be done, using `scalarize` function.

```julia
julia> Symbolics.scalarize(AAt[2,3])
A[2, 1]*A[3, 1] + A[2, 2]*A[3, 2] + A[2, 3]*A[3, 3]

julia> @syms i::Int j::Int
(i, j)

julia> Symbolics.scalarize(AAt[i,j])
A[i, 1]*A[j, 1] + A[i, 2]*A[j, 2] + A[i, 3]*A[j, 3]
```

In general any scalar expression which is derived from array expressions can be scalarized.

```julia
julia> sum(A[:,1]) + sum(A[2,:])
Symbolics._mapreduce(identity, +, A[Colon(), 1], Colon(), (:init => false,)) + Symbolics._mapreduce(identity, +, A[2, Colon()], Colon(), (:init => false,))

julia> Symbolics.scalarize(sum(A[:,1]) + sum(A[2,:]))
A[1, 1] + A[2, 2] + A[2, 3] + A[4, 1] + A[5, 1] + 2A[2, 1] + A[3, 1]

```
