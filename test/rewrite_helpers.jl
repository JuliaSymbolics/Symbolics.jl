
# Fetch packages.
using Symbolics
using Symbolics.RewriteHelpers
using Test 

@variables t X(t) Y(t) Z(t)
@syms a b
D = Differential(t)
my_f(x, y) = x^3 + 2y

# Check replace function.
let
    @test isequal(replace(X + X + X, X =>1), 3)
    @test isequal(replace(X + X + X, Y => 1), 3X)
    @test isequal(replace(X + X + X, X => Y), 3Y)
    @test isequal(replace(X + Y^2 - Z, Y^2 => Z), X)
end

# Test occursin function.
let
    ex1 = 2X^a - log(b + my_f(Y,Y)) - 3
    ex2 = X^(Y^(Z-a)) +log(log(log(b)))
    ex3 = sin(X) + sin(Y) + a*a*a*(1-X)
    ex4 = exp(a)/(pi*a) + D(Y) + D(my_f(1,Z))
    ex5 = a + 5b^2
    
    # Test for variables.
    @test occursin(X, ex1)
    @test occursin(X, ex2)
    @test occursin(X, ex3)
    @test !occursin(X, ex4)
    @test occursin(Y, ex1)
    @test occursin(Y, ex2)
    @test occursin(Y, ex3)
    @test occursin(Y, ex4)
    @test !occursin(Z, ex1)
    @test occursin(Z, ex2)
    @test !occursin(Z, ex3)
    @test occursin(Z, ex4)
    
    # Test for variables.
    @test_broken occursin(a, ex1)
    @test_broken occursin(a, ex2)
    @test_broken occursin(a, ex3)
    @test_broken occursin(a, ex4)
    @test occursin(a, ex5)
    @test_broken occursin(b, ex1)
    @test_broken occursin(b, ex2)
    @test !occursin(b, ex3)
    @test !occursin(b, ex4)
    @test occursin(b, ex5)
    
    # Test for function.
    @test !occursin(is_derivative, ex1)
    @test !occursin(is_derivative, ex2)
    @test !occursin(is_derivative, ex3)
    @test occursin(is_derivative, ex4)
end

# Check filterchildren function.
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

    # Test for variables.
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

