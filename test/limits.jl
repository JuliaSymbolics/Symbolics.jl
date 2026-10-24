using Symbolics, Test

@testset "limits" begin
    @syms x
    @test limit(exp(x+exp(-x))-exp(x), x, Inf) == 1
    @test limit(x^7/exp(x), x, Inf) == 0
    @test limit(log(log(x*exp(x*exp(x))+1))-exp(exp(log(log(x))+1/x)), x, Inf) == 0
    @test limit(2exp(-x)/exp(-x), x, 0) == 2

    @test_throws ArgumentError limit(1/x, x, 0)
    @test limit(1/x, x, 0, :left)[1] == -Inf
    @test limit(1/x, x, 0, :right)[1] == Inf
end
