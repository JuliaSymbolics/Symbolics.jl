# Function Building and Compilation (build_function)

At any time, callable functions can be generated from Symbolics IR by
using `Symbolics.toexpr`. This performs some cleaning to return an
expression without extraneous pieces that commonly matches expressions one
would write in functions like those for differential equation solvers and
optimization libraries. These functions can be automatically parallelize and
specialize on Julia types like static arrays and sparse matrices.

The core compilation process of Symbolics IR is `build_function`.
`build_function` takes an operation or an `AbstractArray` of operations and
generates a compilable version of the model for numerical solvers. The form of
this output is dependent on the `target`. By default, the target outputs
Julia code, but other formats, such as C, Stan, and MATLAB are available.
These can be generated as expressions which can then be evaluated into a callable
function, or the compilers for the respective targets can be invoked to directly
give back the function handle.

## build_function

```@docs
build_function
```

## Target-Specific Definitions

```@docs
Symbolics._build_function(target::Symbolics.JuliaTarget,rhss::AbstractArray,args...;kwargs...)
Symbolics._build_function(target::Symbolics.CTarget,eqs::Array{<:Equation},args...;kwargs...)
Symbolics._build_function(target::Symbolics.StanTarget,eqs::Array{<:Equation}, vs, ps, iv;kwargs...)
Symbolics._build_function(target::Symbolics.MATLABTarget,eqs::Array{<:Equation},args...;kwargs...)
```
