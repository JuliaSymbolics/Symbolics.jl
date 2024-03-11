using Test, Symbolics
using Symbolics: Struct, juliatype, symbolic_getproperty, symbolic_setproperty!

struct Jörgen
    a::Int
    b::Float64
end

S = Struct(Jörgen)
@variables x::S
xa = Symbolics.unwrap(symbolic_getproperty(x, :a))
@test Symbolics.symtype(xa) == Int
@test Symbolics.operation(xa) == getfield
@test isequal(Symbolics.arguments(xa), [Symbolics.unwrap(x), Meta.quot(:a)])
xa = Symbolics.unwrap(symbolic_setproperty!(x, :a, 10))
@test Symbolics.operation(xa) == setfield!
@test isequal(Symbolics.arguments(xa), [Symbolics.unwrap(x), Meta.quot(:a), 10])
@test Symbolics.symtype(xa) == Int

xb = Symbolics.unwrap(symbolic_setproperty!(x, :b, 10))
@test Symbolics.operation(xb) == setfield!
@test isequal(Symbolics.arguments(xb), [Symbolics.unwrap(x), Meta.quot(:b), 10])
@test Symbolics.symtype(xb) == Float64
