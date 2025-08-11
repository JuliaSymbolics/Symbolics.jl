# Metadata

Symbolics.jl provides a metadata system for attaching additional information to symbolic variables. This system allows for extensible annotations that can be used to store default values, source information, and custom user-defined metadata.

## Metadata Types

```@docs
Symbolics.VariableDefaultValue
Symbolics.VariableSource
Symbolics.CallWithMetadata
Symbolics.Unknown
```

## Extending Metadata

You can define custom metadata types for use with the `@variables` macro:

```@docs
Symbolics.option_to_metadata_type
```

## Example

```julia
using Symbolics

# Define a custom metadata type
struct MyCustomMetadata <: Symbolics.AbstractVariableMetadata end

# Register it for use in @variables
Symbolics.option_to_metadata_type(::Val{:my_custom}) = MyCustomMetadata

# Now you can use it
@variables x [my_custom = "some value"]
```