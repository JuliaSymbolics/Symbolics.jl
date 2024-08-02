# Solver

`symbolic_solve` is a function which attempts to solve input equations/expressions symbolically using various methods.
It expects 2-3 arguments.

- expr: Could be a single univar expression in the form of a poly or multiple univar expressions or multiple multivar polys or a transcendental nonlinear function which is solved by isolation, attraction and collection.

- x: Could be a single variable or an array of variables which should be solved

- multiplicities (optional): Should the output be printed `n` times where `n` is the number of occurrence of the root? Say we have `(x+1)^2`, we then have 2 roots `x = -1`, by default the output is `[-1]`, If multiplicities is inputted as true, then the output is `[-1, -1]`.

It can take a single variable, a vector of variables, a single expression, an array of expressions.
The base solver (`symbolic_solve`) has multiple solvers which chooses from depending on the the type of input
(multiple/uni var and multiple/single expression) only after ensuring that the input is valid.

Solve's `solve_univar` uses analytic solutions up to polynomials of degree 4 and factoring as its method
for solving univariate polynomials. Solve's `solve_multipoly` uses GCD on the input polynomials then throws passes the result
to `solve_univar`. Solve's `solve_multivar` uses Groebner basis and a separating form in order to create linear equations in the
input variables and a single high degree equation in the separating variable. Each equation resulting from the basis is then passed
to `solve_univar`. We can see that essentially, `solve_univar` is the building block of `solve`.

if the input is not a valid polynomial and can not be solved by the algorithm above, `symbolic_solve` passes
it to `ia_solve`, which attempts solving by attraction and isolation [^1]. This only works when the input is a single expression
and the user wants the answer in terms of a single variable. Say `log(x) - a == 0` gives us `[e^a]`.


## Available solvers
The `solve` function implements the following backends:

- `solve_univar` (single variable single polynomial)
- `solve_multivar` (multiple variables multiple polynomials)
- `solve_multipoly` (single variable multiple polynomials)
- `ia_solve` (not in the form of a polynomial, uses isolation and attraction in order to reshape the expression in the form of a poly)


## Examples and usage

#### `solve_univar`
```jldoctest
julia> expr = expand((x + 3)*(x^2 + 2x + 1)*(x + 2))
6 + 17x + 17(x^2) + 7(x^3) + x^4

julia> symbolic_solve(expr, x)
3-element Vector{Any}:
 -2//1
 -1//1
 -3//1

julia> symbolic_solve(expr, x, true)
4-element Vector{Any}:
 -2//1
 -1//1
 -1//1
 -3//1
```
```jldoctest
julia> symbolic_solve(x^7 - 1, x)
2-element Vector{Any}:
roots_of((1//1) + x + x^2 + x^3 + x^4 + x^5 + x^6)
 1//1
```
When a polynomial is not solvable, solve outputs a roots_of struct which has 2 variables,
poly and the which was attempted to solve for.

#### `solve_multivar`
```jldoctest
julia> eqs = [x+y^2+z, z*x*y, z+3x+y]
3-element Vector{Num}:
 x + z + y^2
       x*y*z
  3x + y + z

julia> symbolic_solve(eqs, [x,y,z])
3-element Vector{Any}:
 Dict{Num, Any}(z => 0//1, y => 0//1, x => 0//1)
 Dict{Num, Any}(z => 0//1, y => 1//3, x => -1//9)
 Dict{Num, Any}(z => -1//1, y => 1//1, x => 0//1)
```


#### `solve_multipoly`
```jldoctest
julia> symbolic_solve([x-1, x^3 - 1, x^2 - 1, (x-1)^20], x)
1-element Vector{Rational{BigInt}}:
 1
```

#### `ia_solve`
```jldoctest
julia> symbolic_solve(2^(x+1) + 5^(x+3), x)
1-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
 (-log(2) + 3log(5) - log(complex(-1))) / (log(2) - log(5))
```
```jldoctest
julia> symbolic_solve(log(x+1)+log(x-1), x)
2-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
 (1//2)*RootFinding.ssqrt(8.0)
 (-1//2)*RootFinding.ssqrt(8.0)
```
```jldoctest
julia> symbolic_solve(a*x^b + c, x)
((-c)^(1 / b)) / (a^(1 / b))
```
# References
[^1]: [R. W. Hamming, Coding and Information Theory, ScienceDirect, 1980](https://www.sciencedirect.com/science/article/pii/S0747717189800070).

