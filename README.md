# Symbolics.jl

WIP: Moving the CAS out of ModelingToolkit.jl into this package!

- `@variables` (no `@parameters`) and operations on them
- `simplify`, `substitute` (reexported from SymbolicUtils)
- `@register` (not exported)
- `gradient`, `jacobian`, `hessian`
- `sparsejacobian`, `sparsehessian`
- `jacobian_sparsity`, `hessian_sparsity`
- `lu`, `det` on Array{Num}
