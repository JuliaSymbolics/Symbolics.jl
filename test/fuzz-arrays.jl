using SymbolicUtils
using Test

@syms A[1:2, 1:3] B[1:3, 1:1] C::Matrix{Number} x[1:2] y[1:3] z::Vector{Number}

AA = rand(2,3); BB = rand(3,2); xx = rand(2); yy = rand(3)

leaf_nodes = [A, B, C, x, y, z, AA, BB, xx, yy]

function rand_leaf()
    a = rand(leaf_nodes)
    rand(Bool) ?  a' : a
end

function rand_mul_expr(a=rand_leaf(), b=rand_leaf())
    @show a
    @show b

    try
        @show "hi"
        a * b
        @show size(a * b)
        @show "hello"
    catch err1
        @show "sure enough"
        try
            rand(size(a)...) * rand(size(b)...)
        catch err2

            @test typeof(err1) == typeof(err2)
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
        @test size(a * b) == size(rand(size(a)...) * rand(size(b)...))
    end
end
