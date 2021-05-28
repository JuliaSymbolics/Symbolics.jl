using Symbolics
using SymbolicUtils, Test
using Symbolics: symtype, shape, Unknown, Arr
using Base: Slice
using SymbolicUtils: Sym

@testset "arrays" begin
    @variables X[1:5, 1:5] Y[1:5, 1:5]
    @test typeof(X) <: Arr
    @test shape(X) == Slice.((1:5, 1:5))
    @test shape(Y) == Slice.((1:5, 1:5))

    A = Y[2, :] # Maybe it should be Array vv
    @test typeof(A) <: Arr{Real, 2}
    @test axes(A) == (1:5,)

    B = A[3:5]
    @test axes(B) == (Slice(1:3),)

    i = Sym{Int}(:i)
    j = Sym{Int}(:j)
    @test symtype(X[i,j]) == Real
    @test symtype(X[1,j]) == Real
end
