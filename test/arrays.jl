using Symbolics
using SymbolicUtils, Test
using Symbolics: symtype, shape, wrap, unwrap, Unknown, Arr, arrterm
using Base: Slice
using SymbolicUtils: Sym, term, gethead

@testset "arrays" begin
    @variables X[1:5, 1:5] Y[1:5, 1:5]
    @test_throws BoundsError X[1000]
    @test typeof(X) <: Arr
    @test shape(X) == Slice.((1:5, 1:5))
    @test shape(Y) == Slice.((1:5, 1:5))

    A = Y[2, :]
    @test typeof(A) <: Arr{Num, 1}
    @test axes(A) == (1:5,)

    B = A[3:5]
    @test axes(B) == (Slice(1:3),)

    i = Sym{Int}(:i)
    j = Sym{Int}(:j)
    @test symtype(X[i,j]) == Real
    @test symtype(X[1,j]) == Real
end

@testset "getindex" begin
    @variables X[1:5, 1:5] Y[1:5, 1:5]

    @test isequal(X[1,1], wrap(term(getindex, unwrap(X), 1,1)))

    XX = unwrap(X)
    @test isequal(unwrap(X[1, :]), Symbolics.@arrayop(XX[1,:], (j,), XX[1, j]))
    @test isequal(unwrap(X[:, 2]), Symbolics.@arrayop(XX[:,2], (i,), XX[i, 2]))
    @test isequal(unwrap(X[:, 2:3]), Symbolics.@arrayop(XX[:,2:3], (i,j), XX[i, j], (+), (j in 2:3)))

    @variables t x[1:4](t)
    @syms i::Int
    @test isequal(x[i], gethead(unwrap(x[i]))(t))
end

getdef(v) = getmetadata(v, Symbolics.VariableDefaultValue)
@testset "broadcast & scalarize" begin
    @variables A[1:5,1:3]=42 b[1:3]=[2, 3, 5] t x[1:4](t)
    AA = Symbolics.scalarize(A)
    bb = Symbolics.scalarize(b)
    @test all(isequal(42), getdef.(AA))
    @test getdef.(bb) == [2, 3, 5]
    @test isequal(Symbolics.scalarize([b.*1; b.*1]), [bb; bb])
    @test isequal(Symbolics.scalarize(b.^1), bb)
    c = A*b

    @test isequal(collect(sin.(x)),
                  sin.([x[i] for i in 1:4]))

    @test isequal(Symbolics.scalarize(sin.(A.*c)[1,1]),
                  sin(A[1, 1]*(b[1]*A[1, 1] +
                               b[2]*A[1, 2] +
                               b[3]*A[1, 3])))

    D = Differential(t)
    @test isequal(collect(D.(x) ~ x), map(i->D(x[i]) ~ x[i], eachindex(x)))
    @test_throws ArgumentError A ~ t
end
