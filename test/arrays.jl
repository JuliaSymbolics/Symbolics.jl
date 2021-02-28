using SymbolicUtils, Test
using SymbolicUtils: getmetadata, symtype, ArrayShape, ArrayShapeCtx
using Base: Slice

@testset "arrays" begin
    @syms X[1:5, 2:6] (Y::Real)[1:5, 1:5] Z::Matrix{Float64} i::Int j::Int
    @test symtype(X) == Array{Number, 2}
    @test symtype(Y) == Array{Real, 2}
    @test getmetadata(X, ArrayShapeCtx) == ArrayShape(Slice.((1:5, 2:6)))
    @test getmetadata(Y, ArrayShapeCtx) == ArrayShape(Slice.((1:5, 1:5)))

    A = Y[2, :]
    @test symtype(A) == Array{Real, 1}
    @test axes(getmetadata(A, ArrayShapeCtx)) == (1:5,)

    B = A[3:5]
    @test axes(getmetadata(B, ArrayShapeCtx)) == (Slice(1:3),)

    @test_throws ArgumentError getmetadata(Z[1:2, 3:4], ArrayShapeCtx)

    X[i,j]
    X[1,j]
end
