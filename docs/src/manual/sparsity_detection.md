# Structure and Sparsity Detection

Using the tracing system provided by Symbolics.jl expressions, Symbolics.jl
can automatically detect the sparsity patterns of Julia functions in an efficient
way. This functionality is described in more detail in the paper:

```
@article{gowda2019sparsity,
  title={Sparsity Programming: Automated Sparsity-Aware Optimizations in Differentiable Programming},
  author={Gowda, Shashi and Ma, Yingbo and Churavy, Valentin and Edelman, Alan and Rackauckas, Christopher},
  year={2019}
}
```

Please cite this work if the functionality is used.

## Sparsity Detection
### High-level Dispatches
Calculate the sparsity pattern of the Jacobian or Hessian of a function.
```@docs
Symbolics.jacobian_sparsity(::Function,::AbstractArray,::AbstractArray)
Symbolics.hessian_sparsity(::Function,::AbstractArray)
```
An example on how to use this high-level interface in combination with other Julia packages can be found in the
[SparseDiffTools documentation](https://github.com/JuliaDiff/SparseDiffTools.jl#example).
### Symbolic Dispatches
Calculate the sparsity pattern of the Jacobian or Hessian of a symbolic expression.
```@docs
Symbolics.jacobian_sparsity(::AbstractArray,::AbstractArray)
Symbolics.hessian_sparsity(::Any,::AbstractVector)
```

## Structure Detection

```@docs
Symbolics.islinear
Symbolics.isaffine
```
