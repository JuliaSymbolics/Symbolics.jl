# Discontinuities and Inequality Helpers

These interfaces describe how a discontinuous function behaves at a root or
across a discontinuity. They are used by symbolic solvers and by code that
needs to preserve one-sided or approximate inequality semantics.

```@docs
Symbolics.rootfunction
Symbolics.left_continuous_function
Symbolics.right_continuous_function
Symbolics.majorization_function
Symbolics.minorization_function
Symbolics.approximation_function
Symbolics.@register_discontinuity
```
