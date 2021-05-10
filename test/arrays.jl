using SymbolicUtils, Test
using Symbolics: symtype, shape, Unknown
using Base: Slice

@testset "arrays" begin
    @variables X[1:5, 2:6] (Y::Real)[1:5, 1:5] Z::Matrix{Float64} i::Int j::Int
    @test symtype(X) == Array{Number, 2}
    @test symtype(Y) == Array{Real, 2}
    @test shape(X) == Slice.((1:5, 2:6))
    @test shape(Y) == Slice.((1:5, 1:5))

    A = Y[2, :] # Maybe it should be Array vv
    @test symtype(A) == AbstractArray{Real, 1}
    @test axes(A) == (1:5,)

    B = A[3:5]
    @test axes(B) == (Slice(1:3),)

    @test shape(Z[1:2, 3:4]) isa Unknown

    @test symtype(X[i,j]) == Number
    @test symtype(X[1,j]) == Number
end
