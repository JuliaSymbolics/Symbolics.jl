using Symbolics
using Test
using LambertW


#Testing

@testset "solving tests" begin
    function hasFloat(expr)#make sure answer does not contain any strange floats
        if expr isa Float64
            return !isinteger(expr) && expr != float(pi) && expr != exp(1.0)
        elseif expr isa Equation
            return hasFloat(expr.lhs) || hasFloat(expr.rhs)
        elseif iscall(expr)
            elements = arguments(expr)
            for element in elements
                if hasFloat(element)
                    return true
                end
            end
        end
        return false
    end
    correctAns(p,a) = isapprox(sort(Symbolics.convert_solutions_to_floats(p)), a) && !hasFloat(p)

    @syms x y z a b c

    #quadratics
    @test correctAns(solve_single_eq(x^2~4,x),[-2.0,2.0])
    @test correctAns(solve_single_eq(x^2~2,x),[-sqrt(2.0),sqrt(2.0)])
    @test correctAns(solve_single_eq(x^2~32,x),[-sqrt(32.0),sqrt(32.0)])
    @test correctAns(solve_single_eq(x^3~32,x),[32.0^(1.0/3.0)])
    #lambert w
    @test correctAns(solve_single_eq(x^x~2,x),[log(2.0)/lambertw(log(2.0))])
    @test correctAns(solve_single_eq(2*x*exp(x)~3,x),[LambertW.lambertw(3.0/2.0)])
    #more challenging quadratics
    @test correctAns(solve_single_eq(x+sqrt(1+x)~5,x),[3.0])
    @test correctAns(solve_single_eq(2*x^2-6*x-7~0,x),[(3.0/2.0)-sqrt(23.0)/2.0,(3.0/2.0)+sqrt(23.0)/2.0])
    #functions inverses
    @test correctAns(solve_single_eq(exp(x^2)~7,x),[-sqrt(log(7.0)),sqrt(log(7.0))])
    @test correctAns(solve_single_eq(sin(x+3)~1//3,x),[asin(1.0/3.0)-3.0])
    r = solve_single_eq(sin(x+2//5)+cos(x+2//5)~1//2,x)
    if r !== nothing
        @test correctAns(r, [acos(0.5/sqrt(2.0))+3.141592653589793/4.0-(2.0/5.0)])
    end
    #product
    @test correctAns(solve_single_eq((x^2-4)*(x+1)~0,x),[-2.0,-1.0,2.0])
end
