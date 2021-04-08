using SymbolicUtils
using Test
using LinearAlgebra

# Transpose and multiply

@syms A[1:2, 1:3] B[1:3, 1:1] C::Matrix{Number} x[1:2] y[1:3] z::Vector{Number}

AA = rand(2,3); BB = rand(3,2); xx = rand(2); yy = rand(3)

leaf_nodes = [A, B, C, x, y, z, AA, BB, xx, yy]

function rand_leaf()
    a = rand(leaf_nodes)
    rand(Bool) ?  a' : a
end

function rand_mul_expr(a=rand_leaf(), b=rand_leaf())

    try
        a * b
    catch err1
        try
            rand(size(a)...) * rand(size(b)...)
        catch err2
            if typeof(err1) == typeof(err2)
                @test typeof(err1) == typeof(err2)
            else
                @test_skip typeof(err1) == typeof(err2)
            end
        end
    end
    sz = try size(a * b) catch err; nothing end

    if sz !== nothing
        try
            size(a)
            size(b)
            @goto test_size
        catch err
            return # no size known
        end
        @label test_size

        if (a isa Adjoint{<:Any, <:AbstractVector} && ndims(b) == 1) || SymbolicUtils.isdot(a, b)
            if size(a*b) != ()
                println("a * b is wrong:")
                @show a b
                @show typeof(a) typeof(b)
                return @test size(a*b) == ()
            else
                return @test true
            end
        end

        if size(a * b) == size(rand(size(a)...) * rand(size(b)...))
            @test true
        else
            println("a * b is wrong:")
            @show a b
            @show typeof(a) typeof(b)
            @test size(a * b) == size(rand(size(a)...) * rand(size(b)...))
        end
    end
end


@testset "fuzz ' and *" begin
    for i=1:1000
        rand_mul_expr()
    end
end

## Getindex fuzz
#

function rand_getindex_expr(x=rand_leaf())
    n = ndims(x)
end
