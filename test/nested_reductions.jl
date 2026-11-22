using Symbolics, Test, LinearAlgebra

@testset "Nested lazy array reductions" begin
    @variables x[1:2] p[1:2] y
    xs = collect(x)
    ps = collect(p)
    A = [1 2; 3 4]
    @test isequal(Symbolics.gradient(2sum(A * x), xs), Num[8, 12])
    for scale in (1, 2, 16)
        B = scale * A
        residual = B * xs - ps
        obj = sum(abs2, B * x - p) + 3sum(abs2, x)
        expected = 2B' * residual + 6xs
        @test all(iszero, simplify.(Symbolics.gradient(obj, xs) - expected; expand = true))
        @test all(iszero, simplify.(Symbolics.hessian(obj, xs) - (2B' * B + 6I)))
        @test all(iszero, simplify.(Symbolics.gradient(obj, ps) + 2residual; expand = true))
        @test isequal(Symbolics.derivative(obj, y), 0)
    end
end

@testset "Linearity through lazy reductions" begin
    @variables x[1:2] y
    xs = collect(x)
    A = [1 2; 3 4]
    @test Symbolics.islinear(sum(x), xs)
    @test Symbolics.islinear(sum(A * x), xs)
    @test Symbolics.isaffine(sum(A * x) + 1, xs)
    @test !Symbolics.islinear(sum(A * x) + 1, xs)
    @test !Symbolics.isaffine(sum(abs2, x), xs)
    @test !Symbolics.islinear(sum(A * x) * y, [xs; y])
end

@testset "Lazy reductions as coordinates" begin
    @variables x[1:2] y
    s = sum(x)
    @test Matrix(Symbolics.hessian_sparsity(s^2, [s])) == [true;;]
    @test Matrix(Symbolics.hessian_sparsity(s * y, [s, y])) == [false true; true false]
    @test !Symbolics.isaffine(s^2, [s])
    @test Symbolics.islinear(s, [s])
    @test Symbolics.islinear(s + 1, [s + 1])

    vars = [s, y, x[1], x[2]]
    expected = [false true false false; true false false false; false false true false; false false false true]
    @test Matrix(Symbolics.hessian_sparsity(s * y + sum(abs2, x), vars)) == expected
end

function independent_reduction_derivative(n)
    @variables x[1:n] a[1:n] y
    expr = sqrt(sum(abs2, x .- a))
    bytes = @allocated result = Symbolics.derivative(expr, y)
    return result, bytes
end

@testset "Derivative of an independent lazy reduction" begin
    independent_reduction_derivative(4)
    for n in (256, 2048)
        result, bytes = independent_reduction_derivative(n)
        @test iszero(result)
        @test bytes < 64 * 1024
    end
end

@testset "Reductions over arrays of unknown size" begin
    @variables u::Vector{Real} x[1:2] y
    xs = collect(x)
    A = [1 2; 3 4]
    @test isequal(Symbolics.derivative(sum(u)^2 + y^2, y), 2y)
    @test isequal(Symbolics.gradient(sum(abs2, x) + sum(u)^2, xs), 2xs)
    @test Matrix(Symbolics.hessian_sparsity(y^2 + sum(u), [y])) == [true;;]
    @test Matrix(Symbolics.hessian_sparsity(sum(abs2, A * x) + sum(u), xs)) == trues(2, 2)
    @test Symbolics.islinear(sum(u) * y, [y])
    @test !Symbolics.islinear(sum(u) + y, [y])
    @test Symbolics.isaffine(sum(A * x) + sum(u), xs)
    @test !Symbolics.islinear(sum(A * x) + sum(u), xs)
    @test_throws ErrorException Symbolics.hessian_sparsity(dot(u, A * x)^2, xs)
end
