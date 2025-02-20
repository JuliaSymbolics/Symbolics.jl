using Symbolics
import Symbolics: symbolic_to_float, var_from_nested_derivative, unwrap, 
                  isblock, flatten_expr!, build_expr, get_variables, 
                  is_singleton, diff2term, tosymbol, lower_varname, 
                  makesubscripts, degree, coeff

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

    @variables c(..)
    ex3 = c(x) + c(t) - c(c(t) + y)
    vars3 = Symbolics.get_variables(ex3)
    @test length(vars3) == 4

    sorted_vars3 = Symbolics.get_variables(ex3; sort = true)
    @test isequal(sorted_vars3, [c.f, t, x, y])
end

@testset "symbolic_to_float" begin
    @variables x
    @test symbolic_to_float((1//2 * x)/x) isa Rational{Int}
    @test symbolic_to_float((1/2 * x)/x) isa Float64
    @test symbolic_to_float((1//2)*√(279//4)) isa Float64
    @test symbolic_to_float((big(1)//2)*√(279//4)) isa BigFloat
    @test symbolic_to_float((-1//2)*√(279//4)) isa Float64
end

@testset "var_from_nested_derivative" begin
    @variables t x(t) p(..)
    D = Differential(t)
    @test var_from_nested_derivative(x) == (x, 0)
    @test var_from_nested_derivative(D(x)) == (x, 1)
    @test var_from_nested_derivative(p) == (p, 0)
    @test var_from_nested_derivative(D(p(x))) == (p(x), 1)
end

@testset "fixpoint_sub maxiters" begin
    @variables x y
    expr = Symbolics.fixpoint_sub(x, Dict(x => y, y => x))
    @test isequal(expr, x)
    expr = Symbolics.fixpoint_sub(x, Dict(x => y, y => x); maxiters = 9)
    @test isequal(expr, y)
end

@testset "Issue#1342 substitute working on called symbolics" begin
    @variables p(..) x y
    arg = unwrap(substitute(p(x), [p => identity]))
    @test iscall(arg) && operation(arg) == identity && isequal(only(arguments(arg)), x)
    @test unwrap(substitute(p(x), [p => sqrt, x => 4.0])) ≈ 2.0
    arg = Symbolics.fixpoint_sub(p(x), [p => sqrt, x => 2y + 3, y => 1.0 + p(4)])
    @test arg ≈ 3.0
end

# Helper Functions
@variables x y z t u(x, t) v(t) w[1:2]

@testset "isblock" begin
    @test isblock([Expr(:block, :(x + y))]) == true
    @test isblock([Expr(:call, :f, x, y)]) == false
    @test isblock([Expr(:block, :(begin x; y; end))]) == true
    @test isblock([Expr(:return, :(x + y))]) == false
end

@testset "flatten_expr!" begin
    expr = Expr(:block, :(x + y))
    @test flatten_expr!(expr.args) == Any[:(x + y)]

    expr2 = Expr(:block, :(begin x + y; z; end))
    @test flatten_expr!(expr2.args) == Any[:(x + y), :z]
end

@testset "build_expr" begin
    expr = build_expr(:block, [:(x + y), :(y + z)])
    @test expr.head == :block
    @test expr.args == [:(x + y), :(y + z)]
end

@testset "is_singleton" begin
    @test is_singleton(x) == false
    @test is_singleton(sin(x)) == false
    @test is_singleton(u) == false
end

@testset "tosymbol" begin
    expr1 = sin(x)
    @test tosymbol(expr1) == Symbol("sin(x)")

    expr2 = cos(y)
    @test tosymbol(expr2) == Symbol("cos(y)")

    expr4 = u
    @test tosymbol(expr4) == Symbol("u(x, t)")
end

@testset "degree" begin
    expr1 = x^2 + y^3
    @test degree(expr1, x) == 2
    @test degree(expr1, y) == 3

    expr2 = x * y^2 + x^3
    @test degree(expr2, x) == 3
    @test degree(expr2, y) == 2

    expr3 = 1
    @test degree(expr3, x) == 0
end

@testset "coeff" begin
    expr1 = 3x + 2y
    @test coeff(expr1, x) == 3
    @test coeff(expr1, y) == 2
    @test coeff(expr1, x^2) == 0

    expr2 = x^2 + 3x + 2
    @test coeff(expr2, x) == 3
    @test coeff(expr2, x^2) == 1
end

@testset "makesubscripts" begin
    sub1 = makesubscripts(5)
    @test length(sub1) == 5
    @test typeof(sub1[1]) == SymbolicUtils.BasicSymbolic{Int64}

    sub2 = makesubscripts(10)
    @test length(sub2) == 10
end

@testset "diff2term" begin
    @variables x t u(x, t) z(t)
    Dt = Differential(t)
    Dx = Differential(x)

    test_var = x
    result = diff2term(test_var)
    @test result === x

    test_nested_derivative = Dx(Dt(Dt(u)))
    result = diff2term(Symbolics.value(test_nested_derivative))
    @test typeof(result) === Symbolics.BasicSymbolic{Real}
end

@testset "`fast_substitute` inside array symbolics" begin
    @variables x y z
    @register_symbolic foo(a::AbstractArray, b)
    ex = foo([x, y], z)
    ex2 = Symbolics.fixpoint_sub(ex, Dict(y => 1.0, z => 2.0))
    @test isequal(ex2, foo([x, 1.0], 2.0))
end

@testset "`fast_substitute` of subarray symbolics" begin
    @variables p[1:4] q[1:5]
    @test isequal(p[1:2], Symbolics.fast_substitute(p[1:2], Dict()))
    @test isequal(p[1:2], Symbolics.fast_substitute(p[1:2], p => p))
    @test isequal(q[1:2], Symbolics.fast_substitute(p[1:2], Dict(p => q)))
    @test isequal(q[1:2], Symbolics.fast_substitute(p[1:2], p => q))
end

@testset "`fast_substitute` folding `getindex`" begin
    @variables x[1:3]
    @test isequal(Symbolics.fast_substitute(x[1], Dict(unwrap(x) => collect(unwrap(x)))), x[1])
    @test isequal(Symbolics.fast_substitute(x[1], unwrap(x) => collect(unwrap(x))), x[1])
end

@testset "numerator and denominator" begin
    @variables x y
    num_den(x) = (numerator(x), denominator(x))
    @test num_den(x) == (x, 1)
    @test num_den(1/x) == (1, x)
    @test num_den(x/y) == (x, y)
end
