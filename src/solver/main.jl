# This file uses Groebner, so it needs to be added in the future.

"""
    solve(expr, x, multiplicities=false)

Solve equations symbolically.

Solve is a function which equates the input expression/s to 0 and solves for the input variable/s x (variables).
It can take a single variable, a vector of variables, a single expression, an array of expressions.
The base solver has multiple solvers which chooses from depending on the the type of input (multiple/uni var and multiple/single expression)
only after ensuring that the input is valid.

All the examples in the REPL in the documentation of the other solvers can be repeated using RootFinding.solve too.
This makes it more convenient for the user. The solver first checks if the input is a valid polynomial, if not, attempt solving
by attraction and isolation (`ia_solve`) which is inspired by the paper in the References section. This only works when the input is a single expression and the user wants the answer
in terms of a single variable. Say `log(x) - a == 0` gives us `[e^a]` using ia_solve. If the input is anything else in the form of a poly,
the solver uses 3 polynomial solvers appropriately depending on the input. 

# Available solvers

The `solve` function implements the following backends:

- `solve_univar` (single variable single polynomial)
- `solve_multivar` (multiple variables multiple polynomials)
- `solve_multipoly` (single variable multiple polynomials)
- `ia_solve` (not in the form of a polynomial, uses isolation and attraction in order to reshape the expression in the form of a poly)

# Arguments
- expr: Could be a single univar expression in the form of a poly or multiple univar expressions or multiple multivar polys or a transcendental nonlinear function which is solved by isolation, attraction and collection.

- x: Could be a single variable or an array of variables which should be solved

- multiplicities: Should the output be printed `n` times where `n` is the number of occurrence of the root? Say we have `(x+1)^2`, we then have 2 roots `x = -1`, by default the output is `[-1]`, If multiplicities is inputted as true, then the output is `[-1, -1]`.

# Examples

TODO(Alex): add a parametric polynomial solve (perhaps smaller) to show-case `roots_of`:
julia> f = expand((x - a)^2 * (x^10 - 1) * (a*x^2 + b*x * c) * (x^20 + a))
julia> RootFinding.solve(f, x)

## `solve_univar` (uses factoring and analytic solutions up to degree 4)
```jldoctest
julia> solve(x^7 - 1, x)
2-element Vector{Any}:
roots_of((1//1) + x + x^2 + x^3 + x^4 + x^5 + x^6)
 1//1
```
```jldoctest
julia> expr = expand((x + 3)*(x^2 + 2x + 1)*(x + 2))
6 + 17x + 17(x^2) + 7(x^3) + x^4

julia> solve(expr, x)
3-element Vector{Any}:
 -2//1
 -1//1
 -3//1

julia> solve(expr, x, true)
4-element Vector{Any}:
 -2//1
 -1//1
 -1//1
 -3//1
```

## `solve_multivar` (uses Groebner basis and `solve_univar` to find roots)
```jldoctest
julia> eqs = [x+y^2+z, z*x*y, z+3x+y]
3-element Vector{Num}:
 x + z + y^2
       x*y*z
  3x + y + z

julia> solve(eqs, [x,y,z])
3-element Vector{Any}:
 Dict{Num, Any}(z => 0//1, y => 0//1, x => 0//1)
 Dict{Num, Any}(z => 0//1, y => 1//3, x => -1//9)
 Dict{Num, Any}(z => -1//1, y => 1//1, x => 0//1)
```


## `solve_multipoly` (uses GCD between the input polys)
```jldoctest
julia> solve([x-1, x^3 - 1, x^2 - 1, (x-1)^20], x)
1-element Vector{Rational{BigInt}}:
 1
```

## `ia_solve` (solving by isolation and attraction)
```jldoctest
julia> solve(2^(x+1) + 5^(x+3), x)
1-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
 (-log(2) + 3log(5) - log(complex(-1))) / (log(2) - log(5))
```
```jldoctest
julia> solve(log(x+1)+log(x-1), x)
2-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
 (1//2)*RootFinding.ssqrt(8.0)
 (-1//2)*RootFinding.ssqrt(8.0)
```
```jldoctest
julia> solve(a*x^b + c, x)
((-c)^(1 / b)) / (a^(1 / b))
```
# References
[^1]: [R. W. Hamming, Coding and Information Theory, ScienceDirect, 1980](https://www.sciencedirect.com/science/article/pii/S0747717189800070).
"""
function solve(expr, x, multiplicities=false)
    type_x = typeof(x)
    expr_univar = false
    x_univar = false


    if (type_x == Num || type_x == SymbolicUtils.BasicSymbolic{Real})
        x_univar = true
        @assert is_singleton(unwrap(x)) "Expected a variable, got $x"
    else
        for var in x
            @assert is_singleton(unwrap(var)) "Expected a variable, got $x"
        end
    end

    if !(expr isa Vector) 
        expr_univar = true
        expr = expr isa Equation ? expr.lhs - expr.rhs : expr
        check_expr_validity(expr)
    else
        expr = Vector{Any}(expr)
        for i in eachindex(expr)
            expr[i] = expr[i] isa Equation ? expr[i].lhs - expr[i].rhs : expr[i]
            check_expr_validity(expr[i])
            !check_poly_inunivar(expr[i], x) && throw("Solve can not solve this input currently")
        end
        expr = Vector{Num}(expr)
    end


    if x_univar

        sols = []
        if expr_univar
            sols = check_poly_inunivar(expr, x) ? solve_univar(expr, x, multiplicities) : ia_solve(expr, x)
        else
            for i in eachindex(expr)
                !check_poly_inunivar(expr[i], x) && throw("Solve can not solve this input currently")
            end
            sols = solve_multipoly(expr, x, multiplicities)
        end

        sols = map(postprocess_root, sols)
        return sols
    end

    if !expr_univar && !x_univar
        for e in expr
            for var in x
                @assert check_poly_inunivar(e, var) "This system can not be currently solved by solve."
            end
        end

        sols = solve_multivar(expr, x, multiplicities)
        for sol in sols
            for var in x
                sol[var] = postprocess_root(sol[var])
            end
        end

        return sols 
    end
end

"""
    solve_univar(expression, x, mult=false)
This solver uses analytic solutions up to degree 4 to solve univariate polynomials.
It first handles the special case of the expression being of operation `^`. E.g. ```math (x+2)^{20}```.
We solve this by removing the int `20`, then solving the poly ```math x+2``` on its own.
If the parameter mult of the solver is set to true, we then repeat the found roots of ```math x+2```
twenty times before returning the results to the user.

Step 2 is filtering the expression after handling this special case, and then factoring
it using `factor_use_nemo`. We then solve all the factors outputted using the analytic methods
implemented in the function `get_roots` and its children.

# Arguments
- expr: Single symbolics Num or SymbolicUtils.BasicSymbolic expression. This is equated to 0 and then solved. E.g. `expr = x+2`, we solve `x+2 = 0`

- x: Single symbolics variable

- mult: Print repeated roots or not?

# Examples

"""
function solve_univar(expression, x, mult=false)
    args = []
    mult_n = 1
    expression = unwrap(expression)
    expression = expression isa PolyForm ? SymbolicUtils.toterm(expression) : expression

    # handle multiplicities, i.e. (x+1)^20
    if iscall(expression)
        expr = simplify(deepcopy(expression))
        args = arguments(expr)
        operation = SymbolicUtils.operation(expr)
        if isequal(operation, ^) && args[2] isa Int64
            expression = wrap(args[1])
            mult_n = args[2]
        end
    end

    subs, filtered_expr = filter_poly(expression, x)
    coeffs, constant = polynomial_coeffs(filtered_expr, [x])

    degree = sdegree(coeffs, x)
    u, factors = factor_use_nemo(wrap(filtered_expr))
    factors = convert(Vector{Any}, factors)

    factors_subbed = map(factor -> ssubs(factor, subs), factors)
    arr_roots = []

    if degree < 5 && length(factors) == 1
        arr_roots = get_roots(expression, x)

        # multiplicities
        if mult
            og_arr_roots = deepcopy(arr_roots)
            for i = 1:(mult_n-1)
                append!(arr_roots, og_arr_roots)    
            end
        end

        return arr_roots
    end

    if length(factors) != 1
        @assert isequal(expand(filtered_expr - u*expand(prod(factors))), 0)

        for factor in factors_subbed
            roots = solve_univar(factor, x, mult)
            if isequal(typeof(roots), RootsOf)
                push!(arr_roots, roots)
            else
                append!(arr_roots, roots)
            end
        end
    end


    if isequal(arr_roots, [])
        return RootsOf(wrap(expression), wrap(x))
    end

    return arr_roots
end

# You can compute the GCD between a system of polynomials by doing the following:
# Get the GCD between the first two polys,
# and get the GCD between this result and the following index,
# say: solve([x^2 - 1, x - 1, (x-1)^20], x)
# the GCD between the first two terms is obviously x-1,
# now we call gcd_use_nemo() on this term, and the following,
# gcd_use_nemo(x - 1, (x-1)^20), which is again x-1.
# now we just need to solve(x-1, x) to get the common root in this
# system of equations.
function solve_multipoly(polys::Vector, x::Num, mult=false)
    polys = unique(polys)

    if length(polys) < 1
        throw("No expressions entered")
    end
    if length(polys) == 1
        return solve_univar(polys[1], x, mult)
    end

    gcd = gcd_use_nemo(polys[1], polys[2])

    for i in eachindex(polys)[3:end]
        gcd = gcd_use_nemo(gcd, polys[i])
    end
    
    if isequal(gcd, 1)
        @info "Nemo gcd is 1."
        return []
    end

    return solve_univar(gcd, x, mult)
end


function solve_multivar(eqs::Vector, vars::Vector{Num}, mult=false)
    throw("Groebner bases engine is required. Execute `using Groebner` to enable this functionality.")
end
