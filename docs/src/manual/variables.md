# Variable and Equation Types

Symbolics IR mirrors the Julia AST but allows for easy mathematical
manipulation by itself following mathematical semantics. The base of the IR is
the `Sym` type, which defines a symbolic variable. Registered (mathematical)
functions on `Sym`s (or `istree` objects) return an expression that `istree`.
For example, `op1 = x+y` is one symbolic object and `op2 = 2z` is another, and
so `op1*op2` is another tree object. Then, at the top, an `Equation`, normally
written as `op1 ~ op2`, defines the symbolic equality between two operations.

## Types

`Sym`, `Term`, and `FnType` are from [SymbolicUtils.jl](https://juliasymbolics.github.io/SymbolicUtils.jl/api/). Note that in
Symbolics, we always use `Sym{Real}`, `Term{Real}`, and
`FnType{Tuple{Any}, Real}`. To get the arguments of an `istree` object, use
`arguments(t::Term)`, and to get the operation, use `operation(t::Term)`.
However, note that one should never dispatch on `Term` or test `isa Term`.
Instead, one needs to use `SymbolicUtils.istree` to check if `arguments` and
`operation` is defined.

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

By default, the `@variables` macros return Num-wrapped objects to allow
calling functions which are restricted to `Number` or `Real`.

```@example variables
using Symbolics
@variables t x y z(t);
Symbolics.operation(Symbolics.value(x + y))
```

```@example variables
Symbolics.operation(Symbolics.value(z))
```

```@example variables
Symbolics.arguments(Symbolics.value(x + y))
```

Note that Julia converts irrationals — like `π` and `ℯ` — to `Float64`
whenever they are involved in arithmetic with other numbers, including
integers.  An expression like `2π` will be converted to a float immediately,
so an expression like `2π * x` will leave the symbolic `x` multiplied by a
`Float64`.  It may be preferable to have a symbolic representation of `π`
also, which can be achieved with `Num(π)`.  For generic programming, it may be
helpful to simply redefine the variable `π` to be of the same type as some
other argument, as in

```@example variables
function f(x)
    let π = oftype(x, π)
        1 + (2 // 3 + 4π / 5) * x
    end
end
f(t)
```

This will work for any floating-point input, as well as symbolic input.

## Symbolic Control Flow

Control flow can be expressed in Symbolics.jl in the following ways:

  - `IfElse.ifelse(cond,x,y)`: this is a dispatch-able version of the `ifelse`
    function provided by `IfElse.jl` which allows for encoding conditionals in
    the symbolic branches.

## Inspection Functions

```@docs
SymbolicUtils.istree
SymbolicUtils.operation
SymbolicUtils.arguments
```
