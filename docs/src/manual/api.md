# API Reference

This page contains the complete API reference for Symbolics.jl, organized by category. Some functions are documented in their respective manual sections and are linked here for completeness.

## Variable Creation and Manipulation

### Variable Creation
See [Variable and Equation Types](variables.md) for:
- `@variables`
- `Symbolics.variable`
- `Symbolics.variables`

### Variable Utilities
See [Expression Manipulation](expression_manipulation.md) for:
- `Symbolics.get_variables`

Additional utilities:
```@docs
Symbolics.get_variables!
Symbolics.getparent
```

### Variable Parsing
```@docs
Symbolics._parse_vars
```

### Variable Metadata
```@docs
Symbolics.VariableDefaultValue
Symbolics.VariableSource
Symbolics.option_to_metadata_type
```

## Type Wrappers and Utilities

See [Supported types and dispatch in Symbolics](types.md) for:
- `Symbolics.wrap`
- `Symbolics.unwrap`

## Array Operations

See [Symbolic Arrays](arrays.md) for:
- `Symbolics.Arr`
- `Symbolics.scalarize`
- `Symbolics.shape`

## Symbolic Differentiation

See [Derivatives and Differentials](derivatives.md) for:
- `Symbolics.derivative`
- `Symbolics.gradient`
- `Symbolics.jacobian`
- `Symbolics.sparsejacobian`
- `Symbolics.sparsejacobian_vals`
- `Symbolics.hessian`
- `Symbolics.sparsehessian`
- `Symbolics.sparsehessian_vals`

## Internal Types and Constants

```@docs
Symbolics.CallWithMetadata
Symbolics.NAMESPACE_SEPARATOR
Symbolics.Unknown
```

## Equations

See [Variable and Equation Types](variables.md) for:
- `Equation`
- `Base.:~(::Num, ::Num)`

## Expression Manipulation

See [Expression Manipulation](expression_manipulation.md) for documented functions.

## Function Building

See [Build Function](build_function.md) for:
- `Symbolics.build_function`
- `Symbolics._build_function`

## Parsing

See [Parsing](parsing.md) for:
- `Symbolics.parse_expr_to_symbolic`

## Integration

See [Integration](integration.md) for integration functions.

## Solving

See [Solver](solver.md) for solving functions.

## Limits

See [Limits](limits.md) for limit functions.

## Groebner Basis

See [Groebner Basis](groebner.md) for Groebner basis functions.

## Taylor Series

See [Taylor Series](taylor.md) for Taylor series functions.

## Sparsity Detection

See [Sparsity Detection](sparsity_detection.md) for sparsity detection functions.

## SymbolicUtils Re-exports

See [Variable and Equation Types](variables.md) for:
- `SymbolicUtils.iscall`
- `SymbolicUtils.operation`
- `SymbolicUtils.arguments`