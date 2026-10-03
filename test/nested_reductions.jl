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
