
# Fetch packages.
using Symbolics
using Symbolics.RewriteHelpers
using Test 

@variables t X(t) Y(t) Z(t)
@syms a b
D = Differential(t)
my_f(x, y) = x^3 + 2y

# Check `replacenode` function.
let
    # Simple replacements.
    @test isequal(replacenode(X + X + X, X =>1), 3)
    @test isequal(replacenode(X + X + X, Y => 1), 3X)
    @test isequal(replacenode(X + X + my_f(X, Z), X => Y), Y^3 + 2Y + 2Z)
    @test isequal(replacenode(X + Y^2 - Z, Y^2 => Z), X)

    # When the rule is a function.
    rep_func(expr) = Symbolics.is_derivative(expr) ? b : expr
    @test isequal(replacenode(D(X + Y) - log(a*Z), rep_func), b - log(a*Z))
    @test isequal(replacenode(D(Z^2) + my_f(D(X), D(Y)) + Z, rep_func), b^3 + 3b + Z)
    @test isequal(replacenode(X + sin(Y + a) + a, rep_func), X + sin(Y + a) + a)

    # On non-symbolic inputs.
    @test isequal(replacenode(1, X =>2.0), 1)
    @test isequal(replacenode(1, rep_func), 1)
end

# Test `hasnode` function.
let
    ex1 = 2X^a - log(b + my_f(Y,Y)) - 3
    ex2 = X^(Y^(Z-a)) +log(log(log(b)))
    ex3 = sin(X) + sin(Y) + a*a*a*(1-X)
    ex4 = exp(a)/(pi*a) + D(Y) + D(my_f(1,Z))
    ex5 = a + 5b^2
    
    # Test for variables.
    @test hasnode(X, ex1)
    @test hasnode(X, ex2)
    @test hasnode(X, ex3)
    @test !hasnode(X, ex4)
    @test hasnode(Y, ex1)
    @test hasnode(Y, ex2)
    @test hasnode(Y, ex3)
    @test hasnode(Y, ex4)
    @test !hasnode(Z, ex1)
    @test hasnode(Z, ex2)
    @test !hasnode(Z, ex3)
    @test hasnode(Z, ex4)
    
    # Test for variables.
    @test hasnode(a, ex1)
    @test hasnode(a, ex2)
    @test hasnode(a, ex3)
    @test hasnode(a, ex4)
    @test hasnode(a, ex5)
    @test hasnode(b, ex1)
    @test hasnode(b, ex2)
    @test !hasnode(b, ex3)
    @test !hasnode(b, ex4)
    @test hasnode(b, ex5)
    
    # Test for function.
    @test !hasnode(is_derivative, ex1)
    @test !hasnode(is_derivative, ex2)
    @test !hasnode(is_derivative, ex3)
    @test hasnode(is_derivative, ex4)

    # On non symbolic inputs:
    @test !hasnode(X, 1)
    @test !hasnode(a, 1)
    @test !hasnode(is_derivative, 1)
end

# Check `filterchildren` function.
let 
    ex1 = 2X^a - log(b + my_f(Y,Y)) - 3
    ex2 = X^(Y^(Z-a)) +log(log(log(b)))
    ex3 = sin(X) + sin(Y) + a*a*a*(1-X)
    ex4 = exp(a)/(pi*a) + D(Y) + D(my_f(1,Z))
    ex5 = a + 5b^2

    # Test for variables.
    @test isequal(filterchildren(X, ex1), [X])
    @test isequal(filterchildren(X, ex2), [X])
    @test isequal(filterchildren(X, ex3), [X, X])
    @test isequal(filterchildren(X, ex4), [])
    @test isequal(filterchildren(Y, ex1), [Y, Y])
    @test isequal(filterchildren(Y, ex2), [Y])
    @test isequal(filterchildren(Y, ex3), [Y])
    @test isequal(filterchildren(Y, ex4), [Y])
    @test isequal(filterchildren(Z, ex1), [])
    @test isequal(filterchildren(Z, ex2), [Z])
    @test isequal(filterchildren(Z, ex3), [])
    @test isequal(filterchildren(Z, ex4), [Z])

    # Test for syms.
    @test isequal(filterchildren(a, ex1), [a])
    @test isequal(filterchildren(a, ex2), [a])
    @test isequal(filterchildren(a, ex3), [a])
    @test isequal(filterchildren(a, ex4), [a, a])
    @test isequal(filterchildren(a, ex5), [a])
    @test isequal(filterchildren(b, ex1), [b])
    @test isequal(filterchildren(b, ex2), [b])
    @test isequal(filterchildren(b, ex3), [])
    @test isequal(filterchildren(b, ex4), [])
    @test isequal(filterchildren(b, ex5), [b])

    # Test for function.
    @test isequal(filterchildren(is_derivative, ex1), [])
    @test isequal(filterchildren(is_derivative, ex2), [])
    @test isequal(filterchildren(is_derivative, ex3), [])
    @test isequal(filterchildren(is_derivative, ex4), [D(Y), D(my_f(1,Z))])
end

# https://github.com/JuliaSymbolics/Symbolics.jl/issues/1175
let
    @variables w z α::Real β::Real;
    r3 = @rule ~x * +(~~ys) => sum(map(y-> ~x * y, ~~ys));
    @test r3(2 * (w+w+α+β)) isa Num
end