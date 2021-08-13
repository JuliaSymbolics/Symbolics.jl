# Variable and Equation Types

Symbolics IR mirrors the Julia AST but allows for easy mathematical
manipulation by itself following mathematical semantics. The base of the IR is
the `Sym` type, which defines a symbolic variable. Registered (mathematical)
functions on `Sym`s (or `isterm` objects) return an expression that `isterm`.
For example, `op1 = x+y` is one symbolic object and `op2 = 2z` is another, and
so `op1*op2` is another tree object. Then, at the top, an `Equation`, normally
written as `op1 ~ op2`, defines the symbolic equality between two operations.

## Types

`Sym`, `Term`, and `FnType` are from [SymbolicUtils.jl](https://juliasymbolics.github.io/SymbolicUtils.jl/api/). Note that in
Symbolics, we always use `Sym{Real}`, `Term{Real}`, and
`FnType{Tuple{Any}, Real}`. To get the arguments of a `isterm` object use
`getargs(t::Term)`, and to get the operation, use `gethead(t::Term)`.
However, note that one should never dispatch on `Term` or test `isa Term`.
Instead, one needs to use `TermInterface.isterm` to check if `arguments` and
`gethead` is defined.

```@docs
@variables
Equation
Base.:~(::Num, ::Num)
```

## A note about functions restricted to `Number`s

`Sym` and `Term` objects are NOT subtypes of `Number`. Symbolics provides
a simple wrapper type called `Num` which is a subtype of `Real`. `Num` wraps
either a Sym or a Term or any other object, defines the same set of operations
as symbolic expressions and forwards those to the values it wraps. You can use
`Symbolics.value` function to unwrap a `Num`.

By default, the `@variables` macros return Num-wrapped objects so as to allow
calling functions which are restricted to `Number` or `Real`.

```julia
julia> @variables t x y z(t);

julia> Symbolics.gethead(Symbolics.value(x + y))
+ (generic function with 377 methods)

julia> Symbolics.gethead(Symbolics.value(z))
z(::Any)::Real

julia> Symbolics.getargs(Symbolics.value(x + y))
2-element Vector{Sym{Real}}:
 x
 y
```

## Symbolic Control Flow

Control flow can be expressed in Symbolics.jl in the following ways:

- `IfElse.ifelse(cond,x,y)`: this is a dispatch-able version of the `ifelse`
  function provided by `IfElse.jl` which allows for encoding conditionals in
  the symbolic branches.

## Inspection Functions

```@docs
TermInterface.isterm
TermInterface.gethead
TermInterface.arguments
```
