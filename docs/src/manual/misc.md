# Miscellaneous API

This page documents various utility functions and constants that don't fit into other categories.

## Constants

```@docs
Symbolics.NAMESPACE_SEPARATOR
```

## Utility functions

```@docs
Symbolics.linear_expansion
Symbolics.LinearExpander
Symbolics._toexpr_metadata
Symbolics._toexpr_op
```

## Deprecated

These names are still exported so that existing code keeps working, but each one warns on
use and will be removed in the next breaking release.

```@docs
Symbolics.infimum
Symbolics.supremum
```
