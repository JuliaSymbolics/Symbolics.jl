using Symbolics, Groebner, Nemo, Test

@testset "Multivar components" begin
    @variables x y z a
    polynomial = (x - y) * (x - z) * (y - z)
    witnesses = [Dict(x => 1, y => 1, z => 2), Dict(x => 1, y => 2, z => 1), Dict(x => 2, y => 1, z => 1)]
    for expression in (polynomial, expand(polynomial)), vars in ([x, y, z], [z, y, x])
        solutions = symbolic_solve([expression], vars; warns = false)
        @test length(solutions) == 3
        @test all(solution -> issetequal(keys(solution), vars), solutions)
        @test all(solution -> isequal(simplify(expand(substitute(expression, solution))), 0), solutions)
        for witness in witnesses
            @test any(solution -> all(v -> isequal(Symbolics.value(substitute(solution[v], witness)), witness[v]), vars), solutions)
        end
    end

    equations = [x * y, z * (x - 1)]
    solutions = symbolic_solve(equations, [x, y, z]; warns = false)
    @test length(solutions) == 3
    @test all(solution -> all(eq -> isequal(simplify(expand(substitute(eq, solution))), 0), equations), solutions)
    for witness in (Dict(x => 0, y => 2, z => 0), Dict(x => 2, y => 0, z => 0), Dict(x => 1, y => 0, z => 2))
        @test any(solution -> all(v -> isequal(Symbolics.value(substitute(solution[v], witness)), witness[v]), [x, y, z]), solutions)
    end

    solutions = symbolic_solve([a * polynomial], [x, y, z]; warns = false)
    @test length(solutions) == 3
    @test all(solution -> !haskey(solution, a), solutions)
end
