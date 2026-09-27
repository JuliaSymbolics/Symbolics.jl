# [Complex numbers](@id complex_numbers)

Symbolics has two representations of a complex symbolic scalar.

  - [`Symbolics.SymbolicNumber`](@ref) wraps one expression tree whose numeric domain is
    not known to be real. It is atomic: `exp(im * x)` stays a single exponential and a
    variable `z::Complex` stays a single variable.
  - `Complex{Num}` stores the real and imaginary parts as two `Num`s. It is the explicit
    Cartesian form, built with `Complex(re, im)`.

```@example complex
using Symbolics
@variables x::Real y::Real z::Complex
typeof(z), typeof(im * x), typeof(Complex(x, y))
```

## Numeric-domain contract

The wrapper is chosen from the symtype of the underlying expression. A symtype that is a
subtype of `Real` gives a `Num`, and any other `Number` symtype gives a `SymbolicNumber`.
Since `Num <: Real` and `SymbolicNumber <: Number`, the Julia type always states what is
known about the domain.

Operations narrow and widen accordingly. Real-valued functions of a complex argument
return a `Num`:

```@example complex
typeof.((real(z), imag(z), abs(z), angle(z)))
```

Combining a real symbolic value with a complex constant, or substituting a complex value
into a real expression, widens the result to a `SymbolicNumber`:

```@example complex
typeof(x + im), typeof(substitute(x^2, Dict(x => z)))
```

Promotion follows the same rule, so a collection that mixes `Num` and `SymbolicNumber`
has element type `SymbolicNumber`:

```@example complex
[x, z]
```

An equation between complex expressions is a single `Equation`; it is not split into
real and imaginary parts.

## Working with `Complex{Num}`

Arithmetic between `Complex{Num}` values stays in Cartesian form. Mixing a
`Complex{Num}` with a `SymbolicNumber` promotes to `SymbolicNumber`, which then holds the
explicit node `complex(re, im)`. Simplification and expansion interpret that node as
`re + im * im_part`, so the two forms compare equal after expansion:

```@example complex
c = Complex(x, y)
iszero(simplify(Symbolics.SymbolicNumber(c) - (x + im * y); expand = true))
```

To go the other way, split a `SymbolicNumber` into its parts:

```@example complex
Complex(reim(z)...)
```

## Differentiation

Derivatives are computed on the raw expressions, and the result container is chosen from
the results. If every entry is real, the result has element type `Num`; otherwise it has
element type `SymbolicNumber`. This holds for any mix of real and complex expressions and
variables, for dense and sparse Jacobians and Hessians alike:

```@example complex
Symbolics.jacobian([im * x, x * z], [x, z])
```

```@example complex
Symbolics.jacobian([x^2, x * y], [x, y])
```

## API

```@docs
Symbolics.SymbolicNumber
```
