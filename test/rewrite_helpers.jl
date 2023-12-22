# Fetch packages.
using Symbolics
using Symbolics.RewriteHelpers
using Test 

@variables t X(t) Y(t) Z(t) V(t) W(t)
@syms a b c d e f g
D = Differential(t)
my_f(x, y) = x^3 + 2y

# Check replace function.
let
    @test isequal(replace(X + X + X, X =>1), 3)
    @test isequal(replace(X + X + X, Y => 1), 3X)
    @test isequal(replace(X + X + X, X => Y), 3Y)
    @test isequal(replace(X + Y^2 - Z, Y^2 => Z), X)

    @test isequal(replace(a + a + a, a =>1), 3) # Errors.
    @test isequal(replace(a + X + b, a + X =>b), 2b) # Errors.
    @test isequal(replace(log(log(log(a))), a => b + C), log(log(log(b+c)))) # Errors.


end

# Check is_derivative function
let
    # Single expressions.
    @test is_derivative(D) # Errors.
    @test !is_derivative(t) # Errors.
    @test !is_derivative(X) # Errors.
    @test !is_derivative(a) # Errors.
    @test !is_derivative(1) # Errors.

    # Composite expressions.
    @test is_derivative(D(X)) # Errors.
    @test !is_derivative(D(X) + 3) # Errors.
    @test is_derivative(D(X + 2a*Y)) # Errors.
    @test !is_derivative(D(X) + D(Y)) # Errors.
    @test !is_derivative(my_f(X)) # Errors.


end

# Test occursin function.
let
    ex1 = 2X^a - log(b + my_f(Y,Y)) - 3
    ex2 = X^(Y^(Z-a)) +log(log(log(b)))
    ex3 = sin(X) + sin(Y) + a*a*a*(1-X)
    ex4 = exp(a)/(pi*a) + D(Y) + D(my_f(1,Z))
    ex5 = 1

    # Test X
    @test occursin(X, ex1) # Errors
    @test occursin(X, ex2)
    @test occursin(X, ex3) # Errors
    @test !occursin(X, ex4) # Errors
    @test !occursin(X, ex5) # Errors

    # Test Y
    @test occursin(Y, ex1) # Errors
    @test occursin(Y, ex2) # Errors
    @test occursin(Y, ex3) # Errors
    @test occursin(Y, ex4) # Errors
    @test !occursin(Y, ex5) # Errors

    # Test Z
    @test !occursin(Z, ex1) # Errors.
    @test !occursin(Z, ex2) # Errors.
    @test !occursin(Z, ex3) # Errors.
    @test occursin(Z, ex4) # Errors.
    @test !occursin(Z, ex5) # Errors

    # Test a
    @test occursin(a, ex1) # Errors.
    @test occursin(a, ex2) # Errors.
    @test occursin(a, ex3) # Errors.
    @test occursin(a, ex4) # Errors.
    @test !occursin(a, ex5)

    # Test b
    @test occursin(b, ex1) # Errors
    @test occursin(b, ex2) # Errors
    @test !occursin(b, ex3) # Errors
    @test !occursin(b, ex4) # Errors
    @test !occursin(b, ex5)

    # Checks composites.
    @test occursin(a + b, a + b + c) # Should this work?
end

# Check filterchildren function.
let 
    filterchildren(a, a)
    filterchildren(a, a + 1) # Errors.
    filterchildren(X, X + 1) # Errors.
    filterchildren(X, X + a*(a + X)) # Errors.
    filterchildren(X, D(X)) # Errors.
    filterchildren(D, D(X)) # Errors.
    filterchildren(D, D(X+Y)) # Errors.


end

