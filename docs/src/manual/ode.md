# Symbolic ODE Solving

While not all ODEs have an analytical solution, symbolic ODE solving is provided by Symbolics.jl for 
subsets of cases for which known analytical solutions can be obtained. These expressions can then
be merged with other techniques in order to accelerate code or gain a deeper understanding of real-world
systems.

## Merging Symbolic ODEs with Numerical Methods: ModelingToolkit.jl

If you are looking to merge Symbolics.jl manipulations with numerical solvers such as
[DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/), then we highly recommend
checking out the [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/dev/) system. It represents
systems of differentible-algebraic equations (DAEs, and extension to ODEs) for which a sophisticated
symbolic analysis pipeline is used to generate highly efficient code. This `mtkcompile` pipeline
makes use of all tricks from Symbolics.jl and many that are more specific to numerical code generation,
such as Pantelides index reduction and tearing of nonlinear systems, and it will analytically solve
subsets of the ODE system as finds possible. Thus if you are attempting to use Symbolics to pre-solve
some parts of an ODE analytically, we recommend allowing ModelingToolkit.jl to do this optimization.

Note that if ModelingToolkit is able to analytically solve the equation, it will give an `ODEProblem`
where `prob.u0 === nothing`, and then running `solve` on the `ODEProblem` will give a numerical
`ODESolution` object that on-demand uses the analytical solution to generate any plots or other artifacts.
The analytical solution can be investigated symbolically using `observed(sys)`.

## Symbolically Solving ODEs

Currently there is no native symbolic ODE solver. Though there are bindings to SymPy

!!! note
    This area is currently under heavy development. More solvers will be available in the near future.

### SymPy 

```@docs
Symbolics.sympy_ode_solve
```
