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

