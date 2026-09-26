# Metadata

Symbolics.jl provides a metadata system for attaching additional information to symbolic variables. This system allows for extensible annotations that can be used to store default values, source information, and custom user-defined metadata.

## Using Metadata

Metadata can be attached to variables when they are created with the `@variables` macro. Common metadata includes default values and other annotations:

```julia
using Symbolics

# Variable with default value
@variables x=1.0 y=2.0

# Variables with custom metadata (once registered)
@variables z [description="Temperature in Kelvin"]
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

The `domain` option records the values a variable is assumed to take:

```julia
using Symbolics

@variables x [domain = (10, Inf)]      # an interval, as a `(lo, hi)` tuple
@variables y [domain = v -> v > 0]     # a predicate on candidate values
```

Symbolics stores the value as given and does not interpret it, so the key is a place for
downstream packages and user-written rules to read assumptions from. For instance, a rewrite
rule can use it to drop the absolute value that `sqrt(z^2)` would otherwise need:

```julia
using SymbolicUtils

function is_nonnegative(v)
    domain = Symbolics.getmetadata(Symbolics.unwrap(v), Symbolics.VariableDomain, nothing)
    domain isa Tuple && first(domain) >= 0
end

r = @rule sqrt((~z)^2) => is_nonnegative(~z) ? ~z : sqrt((~z)^2)

@variables p [domain = (0, Inf)] q
r(Symbolics.unwrap(sqrt(p^2)))   # p
r(Symbolics.unwrap(sqrt(q^2)))   # sqrt(q^2), left unchanged
```

## Metadata API

```@docs
Symbolics.VariableDefaultValue
Symbolics.VariableSource
Symbolics.VariableDomain
Symbolics.option_to_metadata_type
```

`Symbolics.Unknown` is reexported from
[SymbolicUtils](https://symbolicutils.juliasymbolics.org/api/#SymbolicUtils.Unknown).
