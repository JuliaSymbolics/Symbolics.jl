using Symbolics
using SymbolicUtils, SymbolicUtils.Code
import SymbolicUtils as SU
using Test

struct Record1
    x::Int
    y::String
    z::Vector{Real}
end

@symstruct Record1

@testset "`isequal` and `hash`" begin
    @variables rec1::Record1
    @test isequal(rec1, rec1)
    @test isequal(rec1, SU.unwrap(rec1))
    @test isequal(SU.unwrap(rec1), rec1)
    @test hash(rec1) == hash(SU.unwrap(rec1))
end

@testset "`getproperty` operation and arguments" begin
    @variables rec1::Record1
    ex = SU.unwrap(rec1.x)
    @test operation(ex) === Symbolics.SymbolicGetproperty{Record1, :x}()
    @test isequal(arguments(ex), [SU.unwrap(rec1)])
end

@testset "`promote_symtype`" begin
    @variables rec1::Record1
    @test SU.promote_symtype(operation(SU.unwrap(rec1.x)), Record1) === Int
    @test SU.promote_symtype(operation(SU.unwrap(rec1.y)), Record1) === String
    @test SU.promote_symtype(operation(SU.unwrap(rec1.z)), Record1) === Vector{Real}
end

@testset "Basic record" begin
    @variables rec::Record1

    @test rec isa Symbolics.SymStruct{Record1}
    ex = rec.x
    @test ex isa Num
    @test SU.symtype(SU.unwrap(ex)) === Int

    ex = rec.y
    @test ex isa Symbolics.SymbolicT
    @test SU.symtype(ex) === String

    ex = rec.z
    @test ex isa Symbolics.Arr{Num, 1}
    @test SU.symtype(SU.unwrap(ex)) === Vector{Real}
    @test SU.shape(SU.unwrap(ex)) === SU.Unknown(1)

    @test toexpr(rec.x) == :(rec.x)
    ex = [rec.y, rec.x + rec.z[1]]
    val = eval(quote
        let rec = Record1(1, "abc", [1.0, 2.0])
            $(toexpr(ex))
        end
    end)
    @test val == ["abc", 2.0]
end

struct Record2{T}
    x::Int
    y::String
    z::Vector{T}
end

@symstruct Record2{T} begin
    shape(:z) = [1:3]
end

@testset "Parametric record with specified shape" begin
    @variables rec::Record2{Int}

    ex = rec.z
    @test ex isa Symbolics.Arr{Num, 1}
    @test SU.symtype(SU.unwrap(ex)) === Vector{Int}
    @test SU.shape(SU.unwrap(ex)) == [1:3]

    ex = rec.z[1]
    @test ex isa Num
    @test SU.symtype(SU.unwrap(ex)) === Int
end

@testset "Recursive struct" begin
    @variables rec::Record2{Record2{Record1}}

    ex = rec.z
    @test ex isa Symbolics.Arr{Symbolics.SymStruct{Record2{Record1}}, 1}
    @test SU.symtype(SU.unwrap(ex)) === Vector{Record2{Record1}}
    @test SU.shape(SU.unwrap(ex)) == [1:3]

    ex = rec.z[1]
    @test ex isa Symbolics.SymStruct{Record2{Record1}}
    @test SU.symtype(SU.unwrap(ex)) === Record2{Record1}

    ex = rec.z[1].z
    @test ex isa Symbolics.Arr{Symbolics.SymStruct{Record1}, 1}
    @test SU.symtype(SU.unwrap(ex)) === Vector{Record1}
    @test SU.shape(SU.unwrap(ex)) == [1:3]

    ex = rec.z[1].z[1]
    @test ex isa Symbolics.SymStruct{Record1}
    @test SU.symtype(SU.unwrap(ex)) === Record1

    ex = rec.z[1].z[1].z
    @test ex isa Symbolics.Arr{Num, 1}
    @test SU.symtype(SU.unwrap(ex)) === Vector{Real}
    @test SU.shape(SU.unwrap(ex)) === SU.Unknown(1)

    @test toexpr(rec.z[1].z[2].z[3]) == :($getindex($getindex($getindex(rec.z, 1).z, 2).z, 3))
    @variables rec::Record2{Record1}
    val = eval(quote
        let rec = Record2{Record1}(1, "A",
                                   [Record1(2, "B", [2.0, 3.0]),
                                       Record1(3, "C", [3.0, 4.0]),
                                       Record1(4, "D", [4.0, 5.0])])
            $(toexpr(rec.x + rec.z[1].x + rec.z[2].z[1] + rec.z[3].z[2]))
        end
    end)
    @test val == 1 + 2 + 3.0 + 5.0
end

abstract type AbstractRecord1 end

@symstruct AbstractRecord1 begin
    shape(:x) = [1:3]
end

struct ConcreteRecord1_1 <: AbstractRecord1
    x::Vector{Int}
end

struct ConcreteRecord1_2 <: AbstractRecord1
    x::Vector{Int}
end

@symstruct ConcreteRecord1_2 begin
    shape(:x) = [1:2]
end

record_op1(x::ConcreteRecord1_1) = sum(x.x)
record_op1(x::ConcreteRecord1_2) = prod(x.x)

@register_symbolic record_op1(x::AbstractRecord1)

record_arrop1(x::Vector{<:AbstractRecord1}) = sum(record_op1, x)

@register_symbolic record_arrop1(x::Vector{AbstractRecord1})

@testset "`@symstruct` of simple abstract type" begin
    @test Symbolics.has_symwrapper(AbstractRecord1)
    @test Symbolics.wrapper_type(AbstractRecord1) == SymStruct{<:AbstractRecord1}
    @test Symbolics.wrapper_type(ConcreteRecord1_1) === SymStruct{ConcreteRecord1_1}

    @variables r1::ConcreteRecord1_1 r2::ConcreteRecord1_2
    @test SU.shape(r1.x) == [1:3]
    @test SU.shape(r2.x) == [1:2]
end

@testset "Function registration with symbolic structs" begin
    @variables r1::ConcreteRecord1_1 r2::ConcreteRecord1_2
    ex = record_op1(r1)
    @test operation(SU.unwrap(ex)) === record_op1

    val = eval(quote
        let r1 = ConcreteRecord1_1([2,2,3])
            $(toexpr(ex))
        end
    end)
    @test val == 7

    ex = record_op1(r2)
    @test operation(SU.unwrap(ex)) === record_op1

    val = eval(quote
        let r2 = ConcreteRecord1_2([2,2,3])
            $(toexpr(ex))
        end
    end)
    @test val == 12

    @variables r3[1:2]::ConcreteRecord1_1 r4[1:2]::ConcreteRecord1_2

    ex = record_arrop1(r3)
    @test operation(SU.unwrap(ex)) === record_arrop1
    val = eval(quote
        let r3 = [ConcreteRecord1_1([2,3,4]), ConcreteRecord1_1([3,4,5])]
            $(toexpr(ex))
        end
    end)
    @test val == 21

    ex = record_arrop1(r4)
    @test operation(SU.unwrap(ex)) === record_arrop1
    val = eval(quote
        let r4 = [ConcreteRecord1_2([2,3,4]), ConcreteRecord1_2([3,4,5])]
            $(toexpr(ex))
        end
    end)
    @test val == 84

    ex = record_arrop1([r3[1], r1])
    @test operation(SU.unwrap(ex)) === record_arrop1
    ex = record_arrop1([r4[1], r2])
    @test operation(SU.unwrap(ex)) === record_arrop1
end

abstract type AbstractRecord2{T} end

@symstruct AbstractRecord2

@testset "Registering type without parameters works" begin
    @test Symbolics.has_symwrapper(AbstractRecord2)
    @test Symbolics.wrapper_type(AbstractRecord2) == SymStruct{<:AbstractRecord2}
    @test Symbolics.has_symwrapper(AbstractRecord2{Int})
    @test Symbolics.wrapper_type(AbstractRecord2{Int}) == SymStruct{<:AbstractRecord2{Int}}
end

@testset "`issymstruct`" begin
    @variables x rec1::Record1
    # SymStruct instance → true
    @test Symbolics.issymstruct(rec1)
    # Unwrapped BasicSymbolic whose symtype is a registered struct → true
    @test Symbolics.issymstruct(SU.unwrap(rec1))
    # Plain scalar variable → false
    @test !Symbolics.issymstruct(x)
    @test !Symbolics.issymstruct(SU.unwrap(x))
    # Field-access result (symtype Int, not a struct type) → false
    @test !Symbolics.issymstruct(SU.unwrap(rec1.x))
    # Non-symbolic value → false
    @test !Symbolics.issymstruct(1)
end

@testset "Linear indexing" begin
    # Record1 has z::Vector{Real} with unknown shape → linear indexing not supported
    @test !Symbolics.symstruct_supports_linear_indexing(Record1)
    # Record2{Int} has shape(:z) = [1:3] → supported
    @test Symbolics.symstruct_supports_linear_indexing(Record2{Int})

    @variables rec::Record2{Int}
    # x::Int (1 element) + y::String (1 element) + z::Vector{Int} shape [1:3] (3 elements) = 5
    @test length(rec) == 5

    @test isequal(rec[1], SU.unwrap(rec.x))
    @test isequal(rec[2], SU.unwrap(rec.y))
    @test isequal(rec[3], SU.unwrap(rec.z[1]))
    @test isequal(rec[4], SU.unwrap(rec.z[2]))
    @test isequal(rec[5], SU.unwrap(rec.z[3]))

    @test_throws BoundsError rec[0]
    @test_throws BoundsError rec[6]

    elems = collect(rec)
    @test length(elems) == 5
    @test all(e -> e isa Symbolics.SymbolicT, elems)
    @test isequal(elems[1], rec[1])
    @test isequal(elems[5], rec[5])
end

@testset "Linear indexing with nested SymStruct" begin
    @test Symbolics.symstruct_supports_linear_indexing(Record2{Record2{Int}})
    # x(1) + y(1) + z::Vector{Record2{Int}} [1:3] × 5 each = 2 + 15 = 17
    @test Symbolics.symstruct_length(Record2{Record2{Int}}) == 17

    @variables rec::Record2{Record2{Int}}
    @test length(rec) == 17

    # z starts at index 3; z[1] contributes indices 3..7 (5 elements)
    @test isequal(rec[3], SU.unwrap(rec.z[1].x))     # z[1].x (first scalar field)
    @test isequal(rec[4], SU.unwrap(rec.z[1].y))     # z[1].y (second scalar field)
    @test isequal(rec[7], SU.unwrap(rec.z[1].z[3]))  # z[1].z[3] (last element of z[1])
    # z[2] starts at index 8
    @test isequal(rec[8], SU.unwrap(rec.z[2].x))
end

@testset "`diff2term` with SymStruct field access" begin
    @variables t rec(t)::Record2{Int}
    D = Differential(t)

    # diff2term(D(rec(t).x)) should equal the .x field of diff2term(D(rec(t)))
    x_term = SU.unwrap(rec.x)
    dt_x = Symbolics.diff2term(D(x_term))
    @test SU.symtype(dt_x) === Int

    rec_term = SU.unwrap(rec)
    dt_rec = Symbolics.diff2term(D(rec_term))
    @test isequal(dt_x, SU.unwrap(getproperty(Symbolics.SymStruct{Record2{Int}}(dt_rec), :x)))

    # diff2term(D(rec(t).z[2])) should equal diff2term(D(rec(t))).z[2]
    z2_term = SU.unwrap(rec.z[2])
    dt_z2 = Symbolics.diff2term(D(z2_term))
    @test SU.symtype(dt_z2) === Int
    z_arr = getproperty(Symbolics.SymStruct{Record2{Int}}(dt_rec), :z)
    @test isequal(dt_z2, SU.unwrap(z_arr[2]))
end
