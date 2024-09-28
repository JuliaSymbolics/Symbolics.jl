using Test, Symbolics
using Symbolics: symstruct, juliatype, symbolic_getproperty, symbolic_setproperty!,
                 symbolic_constructor, BasicSymbolic

struct Jörgen
    a::Int
    b::Float64
end

S = symstruct(Jörgen)
@variables x::S
xa = Symbolics.unwrap(symbolic_getproperty(x, :a))
@test Symbolics.symtype(xa) == Int
@test Symbolics.operation(xa) == Symbolics.typed_getfield
@test isequal(Symbolics.arguments(xa), BasicSymbolic[Symbolics.unwrap(x), Val{:a}()])
xa = Symbolics.unwrap(symbolic_setproperty!(x, :a, 10))
@test Symbolics.operation(xa) == setfield!
@test isequal(Symbolics.arguments(xa), BasicSymbolic[Symbolics.unwrap(x), Meta.quot(:a), 10])
@test Symbolics.symtype(xa) == Int

xb = Symbolics.unwrap(symbolic_setproperty!(x, :b, 10))
@test Symbolics.operation(xb) == setfield!
@test isequal(Symbolics.arguments(xb), BasicSymbolic[Symbolics.unwrap(x), Meta.quot(:b), 10])
@test Symbolics.symtype(xb) == Float64

s = Symbolics.symbolic_constructor(S, 1, 1.0)
@test Symbolics.symtype(s) == S
