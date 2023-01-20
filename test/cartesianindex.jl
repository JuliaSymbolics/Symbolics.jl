using Symbolics, Test
using Symbolics: Arr
#using SymbolicUtils: substitute

@testset "Symbolic CartesianIndex" begin
    @syms i::Int j::Int k::Int
    I = CartesianIndex(i, j, k)
    @test isequal(I[1], i)
    @test isequal(I[2], j)
    @test isequal(I[3], k)

    J = CartesianIndex(1, 2, 3) + I
    @test isequal(J[1], 1 + i)
    @test isequal(J[2], 2 + j)
    @test isequal(J[3], 3 + k)

    @test isequal(I + I, CartesianIndex(2i, 2j, 2k))
    @test isequal(I + CartesianIndex(1, 2, 3), CartesianIndex(i+1, j+2, k+3))
    @test isequal(CartesianIndex(1, 2, 3) + I, CartesianIndex(1+i, 2+j, 3+k))

    @test isequal(I - I, CartesianIndex(0, 0, 0))
    @test isequal(I - CartesianIndex(1, 2, 3), CartesianIndex(i-1, j-2, k-3))
    @test isequal(CartesianIndex(1, 2, 3) - I, CartesianIndex(1-i, 2-j, 3-k))

    @test isequal(2I, CartesianIndex(2i, 2j, 2k))

    A = rand(2, 3, 4)

    @test isequal(substitute(Arr(A)[I], Dict(i=>1, j=>2, k=>3)), A[1, 2, 3])

    II = substitute(I, Dict(i=>1, j=>2, k=>3))

    @test A[II] == A[1, 2, 3]
end
