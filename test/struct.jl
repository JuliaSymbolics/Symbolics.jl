using Test, Symbolics
using Symbolics: StructElement, Struct, operation, arguments

handledtypes = [Int8,
    Int16,
    Int32,
    Int64,
    UInt8,
    UInt16,
    UInt32,
    UInt64,
    Float16,
    Float32,
    Float64]
for t in handledtypes
    @test Symbolics.decodetyp(Symbolics.encodetyp(t)) === t
end

@variables t x(t)
a = StructElement(Int8, :a)
b = StructElement(Int, :b)
s = Struct([a, b])
sa = s.a
sb = s.b
@test operation(sa) === getfield
@test arguments(sa) == Any[s, 1]
@test arguments(sa) isa Any
@test operation(sb) === getfield
@test arguments(sb) == Any[s, 2]
@test arguments(sb) isa Any

sa1 = (setproperty!(s, :a, UInt8(1)))
@test operation(sa1) === setfield!
@test arguments(sa1) == Any[s, 1, UInt8(1)]
@test arguments(sa1) isa Any

sb1 = (setproperty!(s, :b, "hi"))
@test operation(sb1) === setfield!
@test arguments(sb1) == Any[s, 2, "hi"]
@test arguments(sb1) isa Any
