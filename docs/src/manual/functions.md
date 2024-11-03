# [Function Registration and Tracing](@id function_registration)

Function registration is the ability to define new nodes in the symbolic
graph. This is useful because symbolic computing is declarative, i.e.
symbolic computations express *what* should be computed, not *how* it
should be computed. However, at some level someone must describe
how a given operation is computed. These are the primitive functions,
and a symbolic expression is made up of primitive functions.

Symbolics.jl comes pre-registered with a large set of standard mathematical
functions, like `*` and `sin` to special functions like `erf`, and even
deeper operations like DataInterpolations.jl's `AbstractInterpolation`.
However, in many cases you may need to add your own function, i.e. you
may want to give an imperative code and use this to define a new symbolic
code. Symbolics.jl calls the declaration of new declarative primitives
from an imperative function definition **registration**. This page
describes both the registration process and its companion process,
tracing, for interacting with code written in Julia.

## Direct Tracing

Because Symbolics expressions respect Julia semantics, one way
to generate symbolic expressions is to simply place Symbolics
variables as inputs into existing Julia code. For example, the following
uses the standard Julia function for the Lorenz equations to generate
the symbolic expression for the Lorenz equations:

```@example registration
using Symbolics
function lorenz(du,u,p,t)
 du[1] = 10.0(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end
@variables t p[1:3] u(t)[1:3]
du = Array{Any}(undef, 3)
lorenz(du,u,p,t)
du
```
Or similarly:

```@example registration
@variables t x(t) y(t) z(t) dx(t) dy(t) dz(t) σ ρ β
du = [dx,dy,dz]
u = [x,y,z]
p = [σ,ρ,β]
lorenz(du,u,p,t)
du
```

Note that what has been done here is that the imperative Julia code
for the function `lorenz` has been transformed into a declarative
symbolic graph. Importantly, the code of `lorenz` is transformed
into an expression consisting only of primitive registered functions,
things like `*` and `-`, which come pre-registered with Symbolics.jl
This then allows for symbolic manipulation of the expressions, 
allowing things like simplification and operation
reordering done on its generated expressions. 

### Utility and Scope of Tracing

This notably describes one limitation of tracing: tracing only 
works if the expression being traced is composed of already 
registered functions. If unregistered functions, such as calls
to C code, are used, then the tracing process will error.

However, we note that symbolic tracing by definition does not 
guarantee that the exact choices. The symbolic expressions
may re-distribute the arithmetic, simplify out expressions, or
do other modifications. Thus if this function is function is
sensitive to numerical details in its calculation, one would not
want to trace the function and thus would instead register it
as a new primitive function.

For the symbolic system to be as powerful in its manipulations
as possible, it is recommended that the registration of functions
be minimized to the simplest possible set, and thus registration
should only be used when necessary. This is because any code within
a registered function is treated as a blackbox imperative code that
cannot be manipulated, thus decreasing the potential for simplifications.

## Registering Functions

The Symbolics graph only allows registered Julia functions within its type.
All other functions are automatically traced down to registered
functions. By default, Symbolics.jl pre-registers the common functions
utilized in [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl)
and pre-defines their derivatives. However, the user can utilize the
[`@register_symbolic`](@ref) macro to add their function to allowed functions
of the computation graph.

Additionally, [`@register_array_symbolic`](@ref) can be used to define array functions.
For size propagation it's required that a computation of how the sizes are computed is
also supplied.

### Defining Derivatives of Registered Functions

In order for symbolic differentiation to work, an overload of `Symbolics.derivative` is
required. The syntax is `derivative(typeof(f), args::NTuple{i,Any}, ::Val{j})` where
`i` is the number of arguments to the function and `j` is which argument is being
differentiated. So for example:

```julia
function derivative(::typeof(min), args::NTuple{2,Any}, ::Val{1})
    x, y = args
    IfElse.ifelse(x < y, one(x), zero(x))
end
```

is the partial derivative of the Julia `min(x,y)` function with respect to `x`.

!!! note
    Downstream symbolic derivative functionality only work if every partial derivative that
    is required in the derivative expression is defined. Note that you only need to define
    the partial derivatives which are actually computed.

### Registration of Array Functions

Similar to scalar functions, array functions can be registered to define new primitives for
functions which either take in or return arrays. This is done by using the `@register_array_symbolic`
macro. It acts similarly to the scalar function registration but requires a calculation of the
input and output sizes. For example, let's assume we wanted to have a function that computes the
solution to `Ax = b`, i.e. a linear solve, using an SVD factorization. In Julia, the code for this
would be `svdsolve(A,b) = svd(A)\b`. We would create this function as follows:

```@example array_registration
using LinearAlgebra, Symbolics

svdsolve(A, b) = svd(A)\b
@register_array_symbolic svdsolve(A::AbstractMatrix, b::AbstractVector) begin
    size = size(b)
    eltype = promote_type(eltype(A), eltype(b))
end
```

Now using the function `svdsolve` with symbolic array variables will be kept lazy:

```@example array_registration
@variables A[1:3, 1:3] b[1:3]
svdsolve(A,b)
```

Note that at this time array derivatives cannot be defined.

## Registration API

```@docs
@register_symbolic
@register_array_symbolic
```

## Direct Registration API (Advanced, Experimental)

!!! warn

    This is a lower level API which is not as stable as the macro APIs.

In some circumstances you may need to use the direct API in order to define
registration on functions or types without using the macro. This is done
by directly defining dispatches on symbolic objects.

A good example of this is DataInterpolations.jl's interpolations object.
On an interpolation by a symbolic variable, we generate the symbolic
function (the `term`) for the interpolation function. This looks like:

```julia
using DataInterpolations, Symbolics, SymbolicUtils
(interp::AbstractInterpolation)(t::Num) = SymbolicUtils.term(interp, unwrap(t))
```

In order for this to work, it is required that we define the `symtype` for the
symbolic type inference. This is done via:

```julia
SymbolicUtils.promote_symtype(t::AbstractInterpolation, args...) = Real
```

Additionally a symbolic name is required:

```julia
Base.nameof(interp::AbstractInterpolation) = :Interpolation
```

The derivative is defined similarly to the macro case:

```julia
function Symbolics.derivative(interp::AbstractInterpolation, args::NTuple{1, Any}, ::Val{1})
    Symbolics.unwrap(derivative(interp, Symbolics.wrap(args[1])))
end
```

## Inverse function registration

Symbolics.jl allows defining and querying the inverses of functions.

```@docs
inverse
left_inverse
right_inverse
@register_inverse
has_inverse
has_left_inverse
has_right_inverse
```

Symbolics.jl implements inverses for standard trigonometric and logarithmic functions,
as well as their variants from `NaNMath`. It also implements inverses of
`ComposedFunction`s.
