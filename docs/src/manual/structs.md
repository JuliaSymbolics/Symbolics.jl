# Symbolic Structs

It is occasionally useful to represent structs symbolically, and retain the ability to perform
`getproperty` on them. For example, let us consider the following struct representing a point in
2D cartesian coordinates.

```@example symstruct
using Symbolics, SymbolicUtils

struct Point2D{T}
    x::T
    y::T
end
```

To allow this struct to be used symbolically, the [`@symstruct`](@ref) macro must be used:

```@example symstruct
@symstruct Point2D{T}
```

Note that this macro _requires_ that all type parameters of the struct be present in the macro use.
With this, we can create symbolic structs like any other variable:

```@repl symstruct
@variables point::Point2D{Real}
typeof(point)
@assert point.x isa Num
@assert point.y isa Num
```

Subsequently, fields of the struct can be accessed symbolically:

```@repl symstruct
point.x^2 + point.y^2
```

## Array fields

What if we have a general point type capable of representing arbitrary dimensional cartesian
coordinates?

```@example symstruct
struct PointN{T, N}
    x::Vector{T}

    function PointN{T, N}(x::Vector{T}) where {T, N}
        @assert length(x) == N
        new{T, N}(x)
    end
end
```

The standard `@symstruct PointN{T, N}` would work. However, Symbolics cannot infer the length
of field `x` simply from the type `PointN{T, N}`. As a result, the field would be inferred as
an unknown length vector, and allow arbitrary indexing such as `point.x[1234]`. To provide this
information, the `@symstruct` macro allows specifying additional metadata:

```@example symstruct
@symstruct PointN{T, N} begin
    shape(:x) = (1:N,)
end
```

The `shape` option allows specifying the shape of array fields. Here, we're saying that the shape
of field `:x` is `(1:N,)` (i.e. it is a 1-indexed vector of length `N`). The expression can use
type information of the registered struct type. The specific details of all options are documeneted
in the [`@symstruct`](@ref) macro. With the above declaration, we can now have safe symbolic `PointN`.

```@repl symstruct
@variables point::PointN{Real, 3}
point.x[1]^2 + point.x[2]^2 + point.x[3]^2
@assert size(point.x) == (3,)
try
    point.x[4]
catch e
    Base.showerror(e)
end
```

## Abstract types

Often, it is useful to specify that all subtypes of a given abstract type should be considered
symbolic structs. The `@symstruct` macro can be used on abstract types to enable this behavior.

```@example symstruct
abstract type AbstractRecord{T} end

struct Record1{T} <: AbstractRecord{T}
    x::Vector{T}
    y::T
end

struct Record2{T} <: AbstractRecord{T}
    x::Vector{T}
    z::T
end

@symstruct AbstractRecord{T} begin
    shape(:x) = (1:3,)
end
```

As evidenced above, even options such as `shape` can be specified for abstract types. These are
applicable to all subtypes.

```@repl symstruct
@variables rec1::Record1{Real} rec2::Record2{Real}
@assert size(rec1.x) == (3,)
@assert size(rec2.x) == (3,)
rec1.x[1] + rec2.x[2] + rec1.y + rec2.z
```

The `@symstruct` macro roughly follows Julia subtyping behavior. This means that more specific
`@symstruct` declarations override less specific ones.

```@example symstruct
struct Record3{T} <: AbstractRecord{T}
    x::Vector{T}
    w::T
end

@symstruct Record3{T} begin
    shape(:x) = (1:4,)
end

struct Record4{T} <: AbstractRecord{T}
    x::Vector{T}
    f::Vector{T}
end

@symstruct Record4{T} begin
    shape(:f) = (1:5,)
end
```

```@repl symstruct
@variables rec3::Record3{Real} rec4::Record4{Real}
@assert size(rec3.x) == (4,) # Different from the one declared by `AbstractRecord`
@assert size(rec4.x) == (3,) # Falls back to the `AbstractRecord` definition
@assert size(rec4.f) == (5,) # Uses the more specific declaration
```

## Arrays of symbolic structs

Arrays of symbolic structs are required to have a concrete `eltype`

```@repl symstruct
@variables pts[1:3]::PointN{Real, 3}
pts[1].x[1] + pts[2].x[2]
```

## Nested structs

Symbolic struct types can be nested inside each other.

```@repl symstruct
@variables recs::Record1{Record1{Real}}
recs.x[1].x[1] + recs.x[2].y.y
```

## Registered functions of symbolic structs

There is some unique behavior when registering functions where arguments are (arrays of)
symbolic structs. The `@symstruct` macro works by defining some interface functions to allow
using the `Symbolics.SymStruct` wrapper type, similar to how `Num` is a wrapper type. All
registration functions behave the same as they did prior to this feature if the `SymStruct` is
unwrapped (using `SymbolicUtils.unwrap`). The difference in behavior arises from how wrapped
`SymStruct` types are handled. In general, the `@register_symbolic` and `@register_array_symbolic`
macros can only use information from the provided type annotations. They cannot know to declare
methods for `SymStruct{ConcreteFoo}` if the function is registered with `AbstractFoo`. For
example, the following syntax:

```julia
@register_symbolic foofn(foo::AbstractFoo)
```

Will only work for a variable declared as `@variables foo::ConcreteFoo` if `@symstruct AbstractFoo`
is declared, regardless of whether `@symstruct ConcreteFoo` is declared or not. This will always
work for `@variables afoo::AbstractFoo`.

```julia
@register_symbolic foofn2(foos::Vector{AbstractFoo})
```

Suppose we have declared `@variables a::ConcreteFoo b::ConcreteFoo cs[1:3]::ConcreteFoo`.
`foofn2([a, b])` will work only if `@symstruct AbstractFoo` is declared, regardless of
whether `@symstruct ConcreteFoo` is declared or not. The same applies for `foofn2(cs)`.

If you are not using symbolic structs, the registration macros behave exactly as they
did prior to this feature.

## API

```@docs
@symstruct
SymStruct
```
