# Symbolic Structs

Occasionally it may be useful to register a struct as a symbolic type. This would allow one
to write expressions that involve accessing fields of structs, e.g. `A.x + A.y - 2`. These
can be useful when building functions whose arguments are more complicated structs, as 
opposed to vectors or numberse. As an example, let's take the following representation of 
[quaternions](https://en.wikipedia.org/wiki/Quaternion):
```julia
using Symbolics
struct Quaternion
    a::Real
    b::Real
    c::Real
    d::Real
end
```

To create the symbolic quaternion type, we use the function `symstruct`:
```@example struct
SymQuaternion = Symbolics.symstruct(Quaternion)
@variables q1::SymQuaternion q2::SymQuaternion
SymQuaternion
```
The type of the symbolic struct is `Struct{T}`, where `T` is the wrapped type.

The `getproperty` and `setproperty!` functions are replaced with `symbolic_getproperty`
and `symbolic_setproperty!`. So if we want to obtain a symbolic term that corresponds
to accessing the constant coefficient `q1.a`, we'd write
```@example struct
q1a = Symbolics.symbolic_getproperty(q1, :a)
```
Note that the return is a symbolic term that calls the `typed_getfield` function.

Suppose we also want to write a function representing the constraint that our quaternion
lives on the unit sphere. We can do this as follows:
```@example struct
q1a = Symbolics.symbolic_getproperty(q1, :a)
q1b = Symbolics.symbolic_getproperty(q1, :b)
q1c = Symbolics.symbolic_getproperty(q1, :c)
q1d = Symbolics.symbolic_getproperty(q1, :d)
norm_1 = Symbolics.build_function(q1a^2 + q1b^2 + q1c^2 + q1d^2 - 1 == 0, q1; expression = false)
```
Each of these `q1a` represent something that will evaluate to `q1.a` when actual numeric
objects are substituted.

We can check that the following quaternion has norm zero:
```@example struct
q = Quaternion(0.5, 0.5, 0.5, 0.5)
norm_1(q)
```

```@docs
Symbolics.symstruct
Symbolics.juliatype
Symbolics.symbolic_getproperty
Symbolics.symbolic_setproperty!
```
