# Metadata

Symbolics.jl provides a metadata system for attaching additional information to symbolic variables. This system allows for extensible annotations that can be used to store default values, source information, and custom user-defined metadata.

## Using Metadata

Metadata can be attached to variables when they are created with the `@variables` macro. Common metadata includes default values and other annotations:

```julia
using Symbolics, Latexify

# Variable with default value
@variables x=1.0 y=2.0

# Variables with custom metadata (once registered)
@variables z [description="Temperature in Kelvin"]

# Per-variable LaTeX rendering (see also the Latexification section in I/O)
@variables w0 [latexwrapper = s -> raw"\omega_{0}"]
```

## Extending Metadata

You can define custom metadata types for use with the `@variables` macro by defining a new metadata type and registering it:

```julia
using Symbolics

# Define a custom metadata type
struct MyCustomMetadata <: Symbolics.AbstractVariableMetadata end

# Register it for use in @variables
Symbolics.option_to_metadata_type(::Val{:my_custom}) = MyCustomMetadata

# Now you can use it
@variables x [my_custom = "some value"]
```

## Variable domains

The `domain` option records the set of values a variable is assumed to take, as a DomainSets
`Domain`:

```julia
using Symbolics, DomainSets

@variables x [domain = HalfLine()]       # x >= 0
@variables y [domain = (10, Inf)]        # shorthand for `Interval(10, Inf)`
@variables n [domain = Integers()]
```

That is the representation the `x ∈ Interval(...)` pairings of `Symbolics.VarDomainPairing`
already use, so a domain means the same thing wherever it appears. A `(lo, hi)` tuple is
converted to an `Interval`, exactly as `∈` converts it. Anything else is rejected: a condition
no existing `Domain` expresses is written as a `Domain` subtype with a `Base.in` method, which
keeps it composable with `UnionDomain` and the rest of DomainSets.

Nothing in Symbolics reads the key yet. It is there so that downstream packages and
user-written rewrite rules can act on a variable's assumptions — for instance to drop the
absolute value that `sqrt(z^2)` otherwise needs:

```julia
using SymbolicUtils

function nonnegative(v)
    d = Symbolics.getmetadata(Symbolics.unwrap(v), Symbolics.VariableDomain, nothing)
    d === nothing && return false
    applicable(DomainSets.infimum, d) || return false
    return DomainSets.infimum(d) >= 0
end

r = @rule sqrt((~z)^2) => nonnegative(~z) ? ~z : sqrt((~z)^2)

@variables p [domain = HalfLine()]
@variables q [domain = Interval(-1, 1)]

r(Symbolics.unwrap(sqrt(p^2)))   # p
r(Symbolics.unwrap(sqrt(q^2)))   # sqrt(q^2), left unchanged
```

!!! note
    `issubset` is not a dependable way to ask whether a domain lies inside the nonnegatives:
    `issubset(Interval(10, Inf), HalfLine())` returns `false`, and it throws for `Integers()`.
    The `infimum` test above is used instead, and it is deliberately conservative — for a
    domain where `infimum` does not apply it declines to decide rather than deciding wrongly.

## Metadata API

```@docs
Symbolics.VariableDefaultValue
Symbolics.VariableSource
Symbolics.VariableDomain
Symbolics.option_to_metadata_type
```

`Symbolics.Unknown` is reexported from
[SymbolicUtils](https://symbolicutils.juliasymbolics.org/api/#SymbolicUtils.Unknown).
