using Symbolics
import Symbolics: ssqrt, slog, scbrt, symbolic_solve, ia_solve, postprocess_root

@testset "ia_solve without Nemo" begin
    @test Base.get_extension(Symbolics, :SymbolicsNemoExt) === nothing
    @variables x
    roots = ia_solve(log(2 + x), x)
    roots = @test_warn ["Nemo", "required"] ia_solve(log(2 + x^2), x)
    @test operation(roots[1]) == Symbolics.RootsOf
end

using Groebner, Nemo
E = Base.MathConstants.e

function sort_roots(roots)
    return sort(roots, lt = (x,y) -> real(x)==real(y) ? imag(x)<imag(y) : real(x)<real(y))
end

function sort_arr(sols, vars)
    for i in eachindex(sols) 
        sols[i] = convert(Dict{Num, Any}, sols[i])
        for (var, root) in sols[i]
            sols[i][var] = eval(Symbolics.toexpr(root))
        end
    end

    extract_root(c) = (real(c), imag(c))
    sort_root(dict) = Tuple(extract_root(dict[var]) for var in vars)

    # Sorting the answers lexicographically
    sols = sort(sols, by=sort_root)
end

function check_equal(arr1, arr2)
    l1 = length(arr1)
    if l1 != length(arr2)
        return false
    end
    for i = 1:l1
        if !isequal(keys(arr1[i]), keys(arr2[i]))
            return false
        end
        if !all(values(arr1[i]) .≈ values(arr2[i]))
            return false
        end
    end
    return true
end

function check_approx(arr1, arr2)
    l1 = length(arr1)
    if l1 != length(arr2)
        return false
    end
    for i = 1:l1
        if !isequal(keys(arr1[i]), keys(arr2[i]))
            return false
        end
        if !all(isapprox.(values(arr1[i]), values(arr2[i]), atol=1e-6))
            return false
        end
    end
    return true
end

@variables x y z a b c d e s

@testset "Solving in terms of a constant var" begin
    eq = ((s^2 + 1)/(s^2 + 2*s + 1)) - ((s^2 + a)/(b*c*s^2 + (b+c)*s + d))
    calcd_roots = sort_arr(Symbolics.solve_interms_ofvar(eq, s), [a,b,c,d])
    solve_roots = sort_arr(symbolic_solve(eq, [a,b,c,d]), [a,b,c,d])
    known_roots = sort_arr([Dict(a=>1, b=>1, c=>1, d=>1)], [a,b,c,d])
    @test check_approx(calcd_roots, known_roots)   
    @test check_approx(solve_roots, known_roots)   

    eq = (a+b)*s^2 - 2s^2 + 2*b*s - 3*s
    calcd_roots = sort_arr(Symbolics.solve_interms_ofvar(eq, s), [a,b])
    solve_roots = sort_arr(symbolic_solve(eq, [a,b]), [a,b])
    known_roots = sort_arr([Dict(a=>1/2, b=>3/2)], [a,b])
    @test check_approx(calcd_roots, known_roots)   
    @test check_approx(solve_roots, known_roots)   

    eq = (a*x^2+b)*s^2 - 2s^2 + 2*b*s - 3*s + 2(x^2)*(s^3) + 10*s^3
    calcd_roots = sort_arr(Symbolics.solve_interms_ofvar(eq, s), [a,b,x])
    solve_roots = sort_arr(symbolic_solve(eq, [a,b,x]), [a,b,x])
    known_roots = sort_arr([Dict(a=>-1/10, b=>3/2, x=>-im*sqrt(5)), Dict(a=>-1/10, b=>3/2, x=>im*sqrt(5))], [a,b,x])
    @test check_approx(calcd_roots, known_roots)   
    @test check_approx(solve_roots, known_roots)   
end

@testset "Invalid input" begin
    @test_throws AssertionError symbolic_solve(x, x^2)
end

@testset "Nice univar cases" begin
    found_roots = symbolic_solve(1/x^2 ~ 1/y^2 - 2/x^3 * (x-y), x)
    known_roots = Symbolics.unwrap.([y, -2y])
    @test isequal(found_roots, known_roots)
end

@testset "Deg 1 univar" begin
    @test isequal(symbolic_solve(x+1, x), [-1])

    @test isequal(symbolic_solve(2x+1, x), [-1/2])

    @test isequal(symbolic_solve(x, x), [0]) 

    @test isequal(symbolic_solve((x+1)^20, x), [-1])

    @test isequal(Symbolics.get_roots_deg1(x + y^3, x), [-y^3])

    expr = x - Symbolics.term(sqrt, 2)
    @test isequal(symbolic_solve(expr, x)[1], Symbolics.term(sqrt, 2))

    expr = x + im
    @test symbolic_solve(expr, x)[1] == -im
end

@testset "Deg 2 univar" begin
    expr = x^2 + 1
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(Symbolics.get_roots_deg2(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr, x))))
    arr_known_roots = sort_roots([-im, im])
    @test all(arr_get_roots .≈ arr_known_roots)
    @test all(arr_solve_roots .≈ arr_known_roots)

    expr = x^2 + 2x + 10
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(Symbolics.get_roots_deg2(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr, x))))
    arr_known_roots = sort_roots([-1 + 3im, -1 - 3im])
    @test all(arr_get_roots .≈ arr_known_roots)
    @test all(arr_solve_roots .≈ arr_known_roots)

    expr = x^2 - 10x + 25
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(Symbolics.get_roots_deg2(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr, x))))
    arr_known_roots = [5,5]
    @test all(arr_get_roots .≈ arr_known_roots)
    @test all(arr_solve_roots .≈ arr_known_roots)
end

@testset "Deg 3 univar" begin
    expr = x^3 - 2x^2 + x - 2 
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(Symbolics.get_roots_deg3(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr, x))))
    arr_known_roots = sort_roots([2, -im, im])
    @test all(isapprox.(arr_get_roots, arr_known_roots, atol=0.0000001))
    @test all(arr_solve_roots .≈ arr_known_roots)

    expr = x^3 + x^2 + x + 1
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(Symbolics.get_roots_deg3(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr, x))))
    arr_known_roots = sort_roots([-1, -im, im])
    @test all(isapprox.(arr_get_roots, arr_known_roots, atol=0.0000001))
    @test all(arr_solve_roots .≈ arr_known_roots)

    expr = x^3 + 10x
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(Symbolics.get_roots_deg3(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr, x))))
    arr_known_roots = sort_roots([0, -sqrt(10)*im, sqrt(10)*im])
    @test all(arr_get_roots .≈ arr_known_roots)
    @test all(arr_solve_roots .≈ arr_known_roots)
end

@testset "Deg 4 univar" begin
    expr = x^4 + 1
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(Symbolics.get_roots_deg4(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr, x))))
    arr_known_roots = sort_roots(eval.([-(complex(-1))^(1/4),(complex(-1))^(1/4), (complex(-1))^(3/4), -(complex(-1))^(3/4)]))
    @test all(arr_get_roots .≈ arr_known_roots)
    @test all(arr_solve_roots .≈ arr_known_roots)

    expr = x^4 - 3x^2 + 2
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(Symbolics.get_roots_deg4(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr, x))))
    arr_known_roots = sort_roots(eval.([-1, 1, sqrt(2), -sqrt(2)]))
    @test all(arr_get_roots .≈ arr_known_roots)
    @test all(arr_solve_roots .≈ arr_known_roots)

    expr = x^4 - x^3 - 2x^2 + 6x - 4
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(Symbolics.get_roots_deg4(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr, x))))
    arr_known_roots = sort_roots(eval.([-2, 1, 1-im, 1+im]))
    @test all(arr_get_roots .≈ arr_known_roots)
    @test all(arr_solve_roots .≈ arr_known_roots)

    expr = 386314 - 412163*x - 357800*(x^2) + 1029179*(x^3) - 111927*(x^4)
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(Symbolics.get_roots_deg4(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr, x))))
    r_novar = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr))))
    arr_known_roots = sort_roots(eval.([0.5840484 + 0.4176250im, 0.5840484 - 0.4176250im,
    8.788773679354421, -0.76177906049]))
    @test all(isapprox.(arr_known_roots, arr_get_roots, atol=0.00001))
    @test all(isapprox.(arr_known_roots, arr_solve_roots, atol=0.00001))
    @test all(isapprox.(arr_known_roots, r_novar, atol=0.00001))
end

@testset "Complex coeffs univar" begin
    expr = x^4 + sqrt(complex(-2//1))
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr, x))))
    r_novar = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr))))
    arr_known_roots = sort_roots(eval.([-Complex(-1)^(3/8)*2^(1/8), Complex(-1)^(3/8)*2^(1/8), 
        Complex(-1)^(7/8)*2^(1/8), -Complex(-1)^(7/8)*2^(1/8)]))
    @test all(arr_get_roots .≈ arr_known_roots)
    @test all(r_novar .≈ arr_known_roots)

    # standby
    # expr = x^3 + sqrt(complex(-2//1))*x + 2
end

@testset "Multipoly solver" begin
    @test isequal(symbolic_solve([x^2 - 1, x + 1], x)[1], -1)
    @test isequal(symbolic_solve([x^2 - a^2, x + a], x)[1], -a)
    @test isequal(symbolic_solve([x^20 - a^20, x + a], x)[1], -a)
end
@testset "Multivar solver" begin
    @variables x y z
    @test isequal(symbolic_solve([x^4 - 1, x - 2], [x]), [])
    
    # TODO: test this properly
    sol = symbolic_solve([x^3 + 1, x*y^3 - 1], [x, y])

    eqs = [x*y + 2x^2, y^2 -1]
    arr_calcd_roots = sort_arr(symbolic_solve(eqs, [x,y]), [x,y])
    r_novar = sort_arr(symbolic_solve(eqs), [x,y])
    arr_known_roots = sort_arr([Dict(x=>-1//2, y=>1), Dict(x=>0, y=>-1), Dict(x=>0, y=>1), Dict(x=>1//2, y=>-1)], [x,y])
    arr_known_roots = sort_arr(arr_known_roots, [x,y])
    @test check_equal(arr_calcd_roots, arr_known_roots)   
    @test check_equal(r_novar, arr_known_roots)   

    eqs = [x-y-z, x+y-z^2, x^2 + y^2 - 1]
    eqs_in_equation_form = [x-y ~ z, x+y-z^2 ~ 0, x^2 + y^2 ~ 1]
    arr_calcd_roots = sort_arr(symbolic_solve(eqs, [x,y,z]), [x,y,z])
    r_novar = sort_arr(symbolic_solve(eqs_in_equation_form), [x,y,z])
    r_eq = sort_arr(symbolic_solve(eqs_in_equation_form, [x,y,z]), [x,y,z])
    arr_known_roots = sort_arr([Dict(x => 0, y=>1, z=>-1), Dict(x=>1, y=>0, z=>1),
        Dict(x=>(1/2)*(-2-sqrt(2)*im), y=>(1/2)*(-2+sqrt(2)*im), z=>-sqrt(2)*im),
        Dict(x=>(1/2)*(-2+sqrt(2)*im), y=>(1/2)*(-2-sqrt(2)*im), z=>sqrt(2)*im)], [x,y,z])
    @test check_approx(arr_calcd_roots, arr_known_roots)   
    @test check_approx(r_novar, arr_known_roots)   
    @test check_approx(r_eq, arr_known_roots)   

    eqs = [x^2, y, z]
    eqs_in_equation_form = [x^2 ~ 0, y ~ 0, z ~ 0]
    arr_calcd_roots = sort_arr(symbolic_solve(eqs, [x,y,z], dropmultiplicity=false), [x,y,z])
    r_novar = sort_arr(symbolic_solve(eqs, dropmultiplicity=false), [x,y,z])
    r_eq_novar = sort_arr(symbolic_solve(eqs_in_equation_form, dropmultiplicity=false), [x,y,z])
    arr_known_roots = sort_arr([Dict(x=>0, y=>0, z=>0), Dict(x=>0, y=>0, z=>0)], [x,y,z])
    @test check_equal(arr_calcd_roots, arr_known_roots)   
    @test check_equal(r_novar, arr_known_roots)   
    @test check_equal(r_eq_novar, arr_known_roots)   

    eqs = [y^2 - 1, x]
    arr_calcd_roots = sort_arr(symbolic_solve(eqs, [x,y]), [x,y])
    arr_known_roots = sort_arr([Dict(y=>1//1, x=>0//1), Dict(y=>-1//1, x=>0//1)], [x,y])
    @test check_equal(arr_calcd_roots, arr_known_roots)   


    eqs = [x^5 + x, y]
    arr_calcd_roots = sort_arr(symbolic_solve(eqs, [x,y]), [x,y])
    arr_known_roots = sort_arr([Dict(x=>0, y=>0), Dict(x=>-(complex(-1))^(1/4), y=>0),
    Dict(x=>(complex(-1))^(1/4), y=>0), Dict(x=>-(complex(-1))^(3/4), y=>0),
    Dict(x=>(complex(-1))^(3/4), y=>0)], [x,y])
    @test check_approx(arr_calcd_roots, arr_known_roots)   

    @test isequal(symbolic_solve([x*y - 1, y], [x,y]), [])
    @test isequal(symbolic_solve([x+y+1, x+y+2], [x,y]), [])

    eqs = [-1 + y + z + x^2,
           -1 + x + z + y^2,
           -1 + x + y + z^2]
    arr_calcd_roots = sort_arr(symbolic_solve(eqs, [x,y,z]), [x,y,z])
    arr_known_roots = sort_arr([Dict(z=>0, y=>1, x=>0),
    Dict(z=>1, y=>0, x=>0),
    Dict(z=>0, y=>0, x=>1),
    Dict(z=>-1+sqrt(2), y=>-1+sqrt(2), x=>-1+sqrt(2)),
    Dict(z=>-1-sqrt(2), y=>-1-sqrt(2), x=>-1-sqrt(2))], [x,y,z])
    @test check_approx(arr_calcd_roots, arr_known_roots)    

    # cyclic 3
    @variables z1 z2 z3
    eqs = [z1 + z2 + z3, z1*z2 + z1*z3 + z2*z3, z1*z2*z3 - 1]
    sol = symbolic_solve(eqs, [z1,z2,z3])
    backward = [Symbolics.substitute(eqs, s) for s in sol]
    @test all(x -> all(isapprox.(eval(Symbolics.toexpr(x)), 0; atol=1e-6)), backward)

    @variables x y
    eqs = [2332//232*x + 2131232*y - 1//343434, x + y + 1]
    sol = symbolic_solve(eqs, [x,y])
    backward = [Symbolics.substitute(eqs, s) for s in sol]
    @test all(x -> all(isapprox.(eval(Symbolics.toexpr(x)), 0; atol=1e-6)), backward)

    # TODO: broken
    # @variables H1 H2
    # eqs = [-288421779135875//1125899906842624*H1^2 + 1963378034549373//562949953421312*H1 - 4, -288421779135875//844424930131968*H1*H2 + 1963378034549373//562949953421312*H2 - 4]
    # symbolic_solve(eqs, [H1, H2])

    @test isnothing(symbolic_solve([x*y - 1, sin(x)], [x, y]))

    @variables x y z
    boot = 10
    for i in 1:boot
        # at most 4 roots by Bézout's theorem
        rand_eq(xs, d) = rand(-10:10) + rand(-10:10)*x + rand(-10:10)*y + rand(-10:10)*x*y + rand(-10:10)*x^2 + rand(-10:10)*y^2
        eqs = [rand_eq([x,y],2), rand_eq([x,y],2)]
        sol = symbolic_solve(eqs, [x,y])
        backward = [Symbolics.substitute(eqs, s) for s in sol]
        @test all(x -> all(isapprox.(eval(Symbolics.toexpr(x)), 0; atol=1e-6)), backward)
    end

    @test isnothing(symbolic_solve([x^2, x*y, y^2], [x,y], warns=false))
end

@testset "Multivar parametric" begin
    @variables x y a
    @test isequal(symbolic_solve([x + a, a - 1], x), [-1])
    @test isequal(symbolic_solve([x - a, y + a], [x, y]), [Dict(y => -a, x => a)])
    @test isequal(symbolic_solve([x*y - a, x*y + x], [x, y]), [Dict(y => -1, x => -a)])
    @test isequal(symbolic_solve([x*y - a, 1 ~ 3], [x, y]), [])

    @test isnothing(symbolic_solve([x*y - 1, sin(x)], [x, y]))

    @test isnothing(symbolic_solve([x*y - a, sin(x)], [x, y]))

    @variables t w u v
    sol = symbolic_solve([t*w - 1 ~ 4, u + v + w ~ 1], [t,w])
    @test isequal(sol, [Dict(t => -5 / (-1 + u + v), w => 1 - u - v)])

    sol = symbolic_solve([x-y, y-z], [x])
    @test isequal(sol, [z])

    sol = symbolic_solve([x-y, y-z], [x, y])
    @test isequal(sol, [Dict(x=>z, y=>z)])

    sol = symbolic_solve([x + y - z, y - z], [x])
    @test isequal(sol, [0])
end

@testset "Factorisation" begin
    f = 10x
    u, factors = Symbolics.factor_use_nemo(f)
    @test isequal(u, 10) && isequal(factors, [x])

    f = Symbolics.wrap(10)
    u, factors = Symbolics.factor_use_nemo(f)
    @test isequal(u, 10) && isempty(factors)

    f = x^2 - 1
    u, factors = Symbolics.factor_use_nemo(f)
    @test isequal(u, 1) && isequal(expand(u*prod(factors) - f), 0)

    f = expand((x + 1//3) * ((x*y)^2 + 2x*y + y^2) * (x - z))
    u, factors = Symbolics.factor_use_nemo(f)
    @test isequal(expand(u*prod(factors) - f), 0)
end


# Post Process roots #
@testset "Post Process roots" begin
    SymbolicUtils.@syms __x
    __symsqrt(x) = SymbolicUtils.term(ssqrt, x)
    term = SymbolicUtils.term
    @test Symbolics.postprocess_root(2 // 1) == 2 && Symbolics.postprocess_root(2 + 0*im) == 2
    @test Symbolics.postprocess_root(__symsqrt(4)) == 2
    @test isequal(Symbolics.postprocess_root(__symsqrt(__x)^2), __x)


    @test isequal(Symbolics.postprocess_root(term(^, 0, __x)), 0)
    @test_broken isequal(Symbolics.postprocess_root(term(/, __x, 0)), Inf)
    @test Symbolics.postprocess_root(term(^, __x, 0) ) == 1
    @test Symbolics.postprocess_root(term(^, Base.MathConstants.e, 0) ) == 1
    @test Symbolics.postprocess_root(term(^, Base.MathConstants.pi, 1) ) == Base.MathConstants.pi
    @test isequal(Symbolics.postprocess_root(term(^, __x, 1) ), __x)

    x = Symbolics.term(sqrt, 2)
    @test isequal(Symbolics.postprocess_root( expand((x + 1)^4) ), 17 + 12x)
    @test isequal(Symbolics.postprocess_root( x^5 ), 4 * x)

    @test isequal(Symbolics.postprocess_root(Symbolics.term(sqrt, 9//4)), 3//2)
    @test isequal(Symbolics.postprocess_root(Symbolics.term(sqrt, -27//8)), im*3//2*Symbolics.term(sqrt, 3//2))
end


# Filter Poly # For future reference, i think @check_polynomial is not enough of a test.
@testset "Filter poly" begin
    @variables x y z c1 c2 
    poly = x*sqrt(complex(-2)) + 2.23324234
    subs, filtered_poly = Symbolics.filter_poly(poly, x)
    @test Symbolics.check_polynomial(filtered_poly)

    poly = x + 2im
    subs, filtered_poly = Symbolics.filter_poly(poly, x)
    @test Symbolics.check_polynomial(filtered_poly)

    poly = im*x + Symbolics.wrap(Symbolics.term(sqrt, 2))
    subs, filtered_poly = Symbolics.filter_poly(poly, x)
    @test Symbolics.check_polynomial(filtered_poly)

    poly = (1/im)*x + 3*y*z
    subs, filtered_poly = Symbolics.filter_poly(poly, x)
    @test Symbolics.check_polynomial(filtered_poly)

    poly = (x+1)*Symbolics.term(log, 3)
    subs, filtered_poly = Symbolics.filter_poly(poly, x)
    @test Symbolics.check_polynomial(filtered_poly)
end


@testset "n func occ (smart occurrence counter)" begin
    @test Symbolics.n_func_occ(x, x) == 1
    @test Symbolics.n_func_occ(log(x), x) == 1
    @test Symbolics.n_func_occ(log(x) + x, x) == 2

    # standby
    # @test Symbolics.(log(y) + x , x) == 1
    @test Symbolics.n_func_occ(log(a*x) + b, x) == 1
    
    @test Symbolics.n_func_occ(log(x + sin((x^2 + x)/log(x))), x) == 3
    @test Symbolics.n_func_occ(x^2 + x + x^3, x) == 1
    @test Symbolics.n_func_occ(log(x)^2 - 17, x) == 1
    @test Symbolics.n_func_occ(2^(x^2 + x) + 5^(x+3), x) == 2

    expr = log( log(x) + log(x) ) + log( log(x) + log(x) ) - 11
    @test Symbolics.n_func_occ(expr, x) == 1

    # log(2) - 3log(5) + x*log(2) - x*log(5)
    expr = expand((1 + x)*Symbolics.term(log, 2) - (3 + x)*Symbolics.term(log, 5))
    @test Symbolics.n_func_occ(expr, x) == 1
end



@testset "Isolate/Attract solve" begin
    lhs = ia_solve(a*x^b + c, x)[1]
    lhs2 = symbolic_solve(a*x^b + c, x)[1]
    rhs = Symbolics.term(^, -c.val/a.val, 1/b.val) 
    @test_broken isequal(lhs, rhs)

    @test isequal(symbolic_solve(2/x, x)[1], Inf)
    @test isequal(symbolic_solve(x^1.5, x)[1], 0)

    lhs = symbolic_solve(log(a*x)-b,x)[1]
    @test isequal(Symbolics.unwrap(Symbolics.ssubs(lhs, Dict(a=>1, b=>1))), 1E)
    
    expr = x + 2
    lhs = eval.(Symbolics.toexpr.(ia_solve(expr, x)))
    lhs_solve = eval.(Symbolics.toexpr.(symbolic_solve(expr, x)))
    rhs = [-2]
    @test lhs[1] ≈ rhs[1]
    @test lhs_solve[1] ≈ rhs[1]

    expr = sqrt(log(cbrt(x^2)))
    lhs = sort_roots(eval.(Symbolics.toexpr.(ia_solve(expr, x))))
    lhs_solve = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr, x))))
    rhs = sort_roots([1, -1])
    @test all(isapprox.(lhs, rhs, atol=1e-6))
    @test all(isapprox.(lhs_solve, rhs, atol=1e-6))

    expr = log(x+1) + log(x-1)
    @test isapprox(eval(Symbolics.toexpr(symbolic_solve(expr, x)[1])), sqrt(2), atol=1e-6)

    expr = 2^(x+1) + 5^(x+3)
    lhs = ComplexF64.(eval.(Symbolics.toexpr.(ia_solve(expr, x))))
    lhs_solve = eval.(Symbolics.toexpr.(symbolic_solve(expr, x)))
    rhs = [(-im*Base.MathConstants.pi - log(2) + 3log(5))/(log(2) - log(5))]
    @test lhs[1] ≈ rhs[1]
    @test lhs_solve[1] ≈ rhs[1]

    expr = 3*2^(x+3) + 2*5^(x+1)
    lhs = eval.(Symbolics.toexpr.(ia_solve(expr, x)))
    lhs_solve = eval.(Symbolics.toexpr.(symbolic_solve(expr, x)))
    rhs = [(-im*Base.MathConstants.pi + log(5) - log(12))/(log(2) - log(5))]
    @test lhs[1] ≈ rhs[1]
    @test lhs_solve[1] ≈ rhs[1]

    expr = exp(2x)*exp(x^2 + 3) + 3
    lhs = sort_roots(eval.(Symbolics.toexpr.(ia_solve(expr, x))))
    lhs_solve = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr, x))))
    rhs = sort_roots([-1 + sqrt(-2 + im*Base.MathConstants.pi + log(3)),
        -1 - sqrt(-2 + im*Base.MathConstants.pi + log(3))])
    @test all(lhs .≈ rhs)
    @test all(lhs_solve .≈ rhs)

    expr = x/5 + 3x^2
    lhs = sort_roots(eval.(Symbolics.toexpr.(ia_solve(expr, x))))
    lhs_solve = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr, x))))
    rhs = sort_roots([0, -1//15])
    @test all(lhs .≈ rhs)
    @test all(lhs_solve .≈ rhs)


    expr = log(x)^2 + log(x) + 1
    lhs = sort_roots(eval.(Symbolics.toexpr.(ia_solve(expr, x))))
    lhs_solve = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr, x))))
    rhs = sort_roots([E^(-(complex(-1))^(1/3)), E^((complex(-1))^(2/3))])
    @test all(lhs .≈ rhs)
    @test all(lhs_solve .≈ rhs)

    expr = 2log(x^2 + 1)^2 + 3log(x^2 + 1) + 1
    lhs = sort_roots(eval.(Symbolics.toexpr.(ia_solve(expr, x))))
    lhs_solve = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr, x))))
    rhs = sort_roots([-im*ssqrt(1 - 1/E), im*ssqrt(1 - 1/E),
        -im*ssqrt(1 - 1/ssqrt(E)), im*ssqrt(1 - 1/ssqrt(E))])
    @test all(lhs .≈ rhs)
    @test all(lhs_solve .≈ rhs)


    # this has c1 subbed, ia_solve still needs to incl c1
    # much like sin and cos solutions (i.e. infinite solutions)
    expr = 9^(x^2 + 1) + 3^(x^2 + 1) + 2
    lhs = sort_roots(eval.(Symbolics.toexpr.(ia_solve(expr, x))))
    lhs_solve = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr, x))))
    rhs = sort_roots([-ssqrt(-1 + slog(1/2*(-1 - im*ssqrt(7)))/slog(3)),
        ssqrt(-1 + slog(1/2*(-1 - im*ssqrt(7)))/slog(3)),
        ssqrt(-1 + slog(1/2*(-1 + im*ssqrt(7)))/slog(3)),
        -ssqrt(-1 + slog(1/2*(-1 + im*ssqrt(7)))/slog(3))])

    @test all(lhs .≈ rhs)
    @test all(lhs_solve .≈ rhs)

    @testset "Keyword arguments" begin
        expr = sec(x ^ 2 + 4x + 4) ^ 3 - 3
        roots = ia_solve(expr, x)
        @test length(roots) == 6 # 2 quadratic roots * 3 roots from cbrt(3)
        @test length(Symbolics.get_variables(roots[1])) == 1
        _n = only(Symbolics.get_variables(roots[1]))
        vals = substitute.(roots, (Dict(_n => 0),))
        @test all(x -> isapprox(norm(sec(x^2 + 4x + 4) ^ 3 - 3), 0.0, atol = 1e-14), vals)

        roots = ia_solve(expr, x; complex_roots = false)
        @test length(roots) == 2
        # the `n` in `θ + n * 2π`
        @test length(Symbolics.get_variables(roots[1])) == 1
        _n = only(Symbolics.get_variables(roots[1]))
        vals = substitute.(roots, (Dict(_n => 0),))
        @test all(x -> isapprox(norm(sec(x^2 + 4x + 4) ^ 3 - 3), 0.0, atol = 1e-14), vals)

        roots = ia_solve(expr, x; complex_roots = false, periodic_roots = false)
        @test length(roots) == 2
        @test length(Symbolics.get_variables(roots[1])) == 0
        @test length(Symbolics.get_variables(roots[2])) == 0
        vals = eval.(Symbolics.toexpr.(roots))
        @test all(x -> isapprox(norm(sec(x^2 + 4x + 4) ^ 3 - 3), 0.0, atol = 1e-14), vals)
    end
end

@testset "Sqrt case poly" begin
    # f(x) + sqrt(g(x)) + c
    expr = x + sqrt(x+1) - 5
    lhs_ia = ia_solve(expr, x)[1]
    lhs_att = Symbolics.attract_and_solve_sqrtpoly(expr, x)[1]
    lhs_solve = symbolic_solve(expr, x)[1]
    @test all(isequal(answer, 3) for answer in [lhs_ia, lhs_att, lhs_solve])

    expr = x^2 + x + sqrt(x) + 2
    lhs = sort_roots(eval.(Symbolics.toexpr.(ia_solve(expr, x))))
    lhs_solve = sort_roots(eval.(Symbolics.toexpr.(symbolic_solve(expr, x))))
    rhs = sort_roots([-0.860929555 - 1.604034315im, -0.860929555 + 1.604034315im])
    @test all(isapprox.(lhs, rhs, atol=1e-6))
    @test all(isapprox.(lhs_solve, rhs, atol=1e-6))
    @test all(isequal.(lhs, lhs_solve))
end

@testset "Turn to poly" begin
    @variables x
    # does not sub because these can not be solved as polys
    expr = log(x+1)^2 + log(x) + 1
    expr, sub = Symbolics.turn_to_poly(expr, x)
    @test isequal(sub, Dict())

    expr = 2log(x)^2 + sin(x) + 1
    expr, sub = Symbolics.turn_to_poly(expr, x)
    @test isequal(sub, Dict())

    
    # subs and turns to polys
    expr = log(x)^2 + log(x) + 1
    expr, sub = Symbolics.turn_to_poly(expr, x)
    @test Symbolics.check_polynomial(expr)

    expr = 2log(x^2 + 1)^2 + 3log(x^2 + 1) + 1
    expr, sub = Symbolics.turn_to_poly(expr, x)
    @test Symbolics.check_polynomial(expr)

    expr = sin(x^2)^2 + 3sin(x^2) + 1
    expr, sub = Symbolics.turn_to_poly(expr, x)
    @test Symbolics.check_polynomial(expr)

    expr = 9^(x^2 + 1) + 3^(x^2 + 1) + 2
    expr, sub = Symbolics.turn_to_poly(expr, x)
    @test Symbolics.check_polynomial(expr)
end


using LambertW


#Testing

@testset "Algebraic solver tests" begin
    function correctAns(solve_roots, known_roots)
        solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve_roots)))
        known_roots = sort_roots(known_roots)
        r = isapprox.(solve_roots, known_roots, atol=1e-6)
        isequal(r, []) && return false
        return all(r)
    end

    #quadratics
    @test correctAns(symbolic_solve(x^2~4,x), [-2, 2])
    @test correctAns(symbolic_solve(x^2~2,x),[-sqrt(2.0),sqrt(2.0)])
    @test correctAns(symbolic_solve(x^2~32,x),[-sqrt(32.0),sqrt(32.0)])
    @test correctAns([symbolic_solve(x^3~32,x)[1]],[32.0^(1.0/3.0)])

    #lambert w, currently unsupported
    @test_broken correctAns(symbolic_solve(x^x~2, x, warns=false),[log(2.0)/LambertW.lambertw(log(2.0))])
    @test_broken correctAns(symbolic_solve(2*x*exp(x)~3, x, warns=false),[LambertW.lambertw(3.0/2.0)])

    #more challenging quadratics
    @test correctAns(symbolic_solve(x+sqrt(1+x)~5,x),[3.0])
    @test correctAns(symbolic_solve(2*x^2-6*x-7~0,x),[(3.0/2.0)-sqrt(23.0)/2.0,(3.0/2.0)+sqrt(23.0)/2.0])

    #functions inverses
    @test correctAns(symbolic_solve(exp(x^2)~7,x),[-sqrt(log(7.0)),sqrt(log(7.0))])

    root = symbolic_solve(sin(x+3)~1//3,x)[1]
    var = Symbolics.get_variables(root)[1]
    root = Symbolics.ssubs(root, Dict(var=>0))
    @test correctAns([root],[asin(1.0/3.0)-3.0])

    @test_broken correctAns(symbolic_solve(sin(x+2//5)+cos(x+2//5)~1//2 , x, warns=false), [acos(0.5/sqrt(2.0))+pi/4.0-(2.0/5.0)])

    #product
    @test correctAns(symbolic_solve((x^2-4)*(x+1)~0,x),[-2.0,-1.0,2.0])
end

@testset "`NaN` handling in `$fn`" for fn in [ssqrt, slog, scbrt]
    x = fn(NaN)
    @test x isa Real && isnan(x)
    x = fn(NaN + NaN * im)
    @test x isa Complex && isnan(x)
end
