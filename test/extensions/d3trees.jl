using Symbolics
using D3Trees
using Test

@variables x y z

@testset "D3Trees Extension" begin
    # Test simple expression
    expr1 = x + y
    tree1 = D3Tree(expr1)
    @test tree1 isa D3Tree
    
    # Test complex expression
    expr2 = sin(x) * cos(y) + z^2
    tree2 = D3Tree(expr2, init_expand = 5)
    @test tree2 isa D3Tree
    
    # Test nested expression
    expr3 = (x + y)^3 * exp(z)
    tree3 = D3Tree(expr3)
    @test tree3 isa D3Tree
    
    # Test with custom kwargs
    expr4 = x^2 + 2*x*y + y^2
    tree4 = D3Tree(expr4, init_expand = 10)
    @test tree4 isa D3Tree
end