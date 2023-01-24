using Symbolics
using Test
using LinearAlgebra

# Transpose and multiply

@variables A[1:2, 1:3] B[1:3, 1:1] x[1:2] y[1:3] C::(Matrix{Real}) z::(Vector{Complex})

stype(T::Type{<:Number}) = Float64
stype(T::Type{<:Complex}) = Complex{Float64}

AA = rand(2, 3);
BB = rand(3, 2);
xx = rand(2);
yy = rand(3);

adjmul_leaf_nodes = (A, B, C, x, y, z, AA, BB, xx, yy)

function adjmul_rand_leaf()
    a = rand(adjmul_leaf_nodes)
    rand(Bool) ? a' : a
end
function isadjvec(x::Symbolics.ArrayOp)
    x.term.f == adjoint && ndims(Symbolics.arguments(x.term)[1]) == 1
end
isadjvec(x) = false

function rand_mul_expr(a = adjmul_rand_leaf(),
                       b = adjmul_rand_leaf())
    try
        a * b
    catch err1
        try
            rand(stype(eltype(a)), size(a)...) * rand(stype(eltype(b)), size(b)...)
        catch err2
            if typeof(err1) == typeof(err2)
                @test typeof(err1) == typeof(err2)
            else
                @test_skip typeof(err1) == typeof(err2)
            end
        end
    end
    sz = try
        size(a * b)
    catch err
        nothing
    end

    if sz !== nothing
        try
            size(a)
            size(b)
            @goto test_size
        catch err
            return # no size known
        end
        @label test_size

        if (isadjvec(Symbolics.unwrap(a)) && ndims(b) == 1) || Symbolics.isdot(a, b)
            if size(a * b) != ()
                println("a * b is wrong:")
                @show a b
                @show typeof(a) typeof(b)
                return @test size(a * b) == ()
            else
                return @test true
            end
        end

        ab_sample = rand(stype(eltype(a)), size(a)...) * rand(stype(eltype(b)), size(b)...)
        if size(a * b) == size(ab_sample)
            @test true
            @test (eltype(a * b) <: Real && eltype(ab_sample) <: Real) ||
                  (eltype(ab_sample) <: Complex && eltype(a * b) <: Complex)
        else
            println("a * b is wrong:")
            @show a b
            @show typeof(a) typeof(b)
            @test size(a * b) == size(rand(size(a)...) * rand(size(b)...))
        end
    end
end

@testset "fuzz ' and *" begin for i in 1:1000
    rand_mul_expr()
end end

## Mapreduce
#

AA = rand(2, 3);
BB = rand(3, 2);
xx = rand(2);
yy = rand(3);
CC = rand(2, 3, 4)

mapreduce_leaf_nodes = [A, B, x, y, AA, BB, CC, xx, yy]

function mapreduce_rand_leaf()
    a = rand(adjmul_leaf_nodes)
    rand(Bool) ? a' : a
end

function rand_map_reduce(x = mapreduce_rand_leaf())
    n = ndims(x)

    reducedims = rand(1:(n + 1), rand(0:n))

    init = rand([nothing, rand(stype(eltype(x)))])
    function test_run(x)
        if init == nothing
            sum(x, dims = reducedims)
        else
            sum(x, dims = reducedims, init = init)
        end
    end

    @test size(test_run(rand(stype(eltype(x)), size(x)...))) == size(test_run(x))
end
