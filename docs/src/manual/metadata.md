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

## Metadata API

```@docs
Symbolics.VariableDefaultValue
Symbolics.VariableSource
Symbolics.CallWithMetadata
Symbolics.Unknown
Symbolics.option_to_metadata_type
```