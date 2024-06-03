# Structure and Sparsity Detection

Using the tracing system provided by Symbolics.jl expressions, Symbolics.jl
can automatically detect the sparsity patterns of Julia functions efficiently.
This functionality is described in more detail in the paper:

```
@article{gowda2019sparsity,
  title={Sparsity Programming: Automated Sparsity-Aware Optimizations in Differentiable Programming},
  author={Gowda, Shashi and Ma, Yingbo and Churavy, Valentin and Edelman, Alan and Rackauckas, Christopher},
  year={2019}
}
```

Please cite this work if the functionality is used.

## Sparsity Detection

```@docs
Symbolics.jacobian_sparsity
Symbolics.hessian_sparsity
```

## Structure Detection

```@docs
Symbolics.islinear
Symbolics.isaffine
```

## ADTypes.jl interface

```@docs
Symbolics.SymbolicsSparsityDetector
```
