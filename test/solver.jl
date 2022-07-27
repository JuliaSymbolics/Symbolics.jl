using Symbolics
using Test


#Testing

function print_code(problem::String)
	expr = Meta.parse(problem)
	out = eval(expr)
	println("$problem -> $out")
	return out
end

@testset "solving tests" begin
	eval(Meta.parse("@syms x"))
	problems = [
		("solve_single_eq(x^2~4,x)",[-2.0,2.0]),
		("solve_single_eq(x^2~2,x)",[-sqrt(2.0),sqrt(2.0)]),
		("solve_single_eq(x^2~32,x)",[-sqrt(32.0),sqrt(32.0)]),
		("solve_single_eq(x^3~32,x)",[32.0^(1.0/3.0)]),
		("solve_single_eq(x^x~2,x)",[log(2.0)/lambertw(log(2.0))]),
		("solve_single_eq(x+sqrt(1+x)~5,x)",[3.0]),
	]
	for problem in problems
		solutions = sort(Symbolics.convert_solutions_to_floats(print_code(problem[1])))
		@test isequal(solutions,problem[2])
	end
end


