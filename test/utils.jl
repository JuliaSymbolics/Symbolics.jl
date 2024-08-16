using Symbolics

@testset "get_variables" begin
    @variables t x y z(t)

    ex1 = x + y + sin(z)
    vars1 = Symbolics.get_variables(ex1)
    @test length(vars1) == 3
    @test allunique(vars1)

    sorted_vars1 = Symbolics.get_variables(ex1; sort = true)
    @test isequal(sorted_vars1, [x, y, z])

    ex2 = x - y
    vars2 = Symbolics.get_variables(ex2)
    @test length(vars2) == 2
    @test allunique(vars2)

    sorted_vars2 = Symbolics.get_variables(ex2; sort = true)
    @test isequal(sorted_vars2, [x, y])
end
