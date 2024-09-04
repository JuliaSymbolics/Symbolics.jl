Base.:^(a::Complex{<:Real}, b::Num) = Symbolics.Pow(a, b)
"""
    symbolic_solve(expr, x; dropmultiplicity=true, warns=true)

`symbolic_solve` is a function which attempts to solve input equations/expressions symbolically using various methods.

## Arguments
- expr: Could be a single univar expression in the form of a poly or multiple univar expressions or multiple multivar polys or a transcendental nonlinear function.

- x: Could be a single variable or an array of variables which should be solved

- dropmultiplicity (optional): Should the output be printed `n` times where `n` is the number of occurrence of the root? Say we have `(x+1)^2`, we then have 2 roots `x = -1`, by default the output is `[-1]`, If dropmultiplicity is inputted as false, then the output is `[-1, -1]`.

- warns (optional): When invalid expressions or cases are inputted, should the solver warn you of such cases before returning nothing? if this is set to false, the solver returns nothing. By default, warns are set to true.

## Supported input
The base solver (`symbolic_solve`) has multiple solvers which chooses from depending on the the type of input
(multiple/uni var and multiple/single expression) only after ensuring that the input is valid.

The expressions inputted can contain parameters, which are assumed to be transcendental.
A parameter "a" is transcendental if there exists no polynomial P with
rational coefficients such that P(a) = 0. Check the examples section.

Currently, `symbolic_solve` supports
- Linear and polynomial equations (with parameters)
- Systems of linear and polynomials equations (without extra parameters, for now)
- Equations with transcendental functions (with parameters)

## Examples

### `solve_univar` (uses factoring and analytic solutions up to degree 4)
!!! note
    The package `Nemo` is needed in order to use `solve_univar` as well as `solve_multipoly`,
    so executing `using Nemo` as you will see in the following examples is necessary; otherwise,
    the function will throw an error.
```jldoctest
julia> using Symbolics, Nemo;

julia> @variables x a b;

julia> expr = expand((x + b)*(x^2 + 2x + 1)*(x^2 - a))
-a*b - a*x - 2a*b*x - 2a*(x^2) + b*(x^2) + x^3 - a*b*(x^2) - a*(x^3) + 2b*(x^3) + 2(x^4) + b*(x^4) + x^5

julia> symbolic_solve(expr, x)
4-element Vector{Any}:
 -1
   -b
   (1//2)*√(4a)
   (-1//2)*√(4a)

julia> symbolic_solve(expr, x, dropmultiplicity=false)
5-element Vector{Any}:
 -1
 -1
   -b
   (1//2)*√(4a)
   (-1//2)*√(4a)
```
```jldoctest
julia> symbolic_solve(x^2 + a*x + 6, x)
2-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
 (1//2)*(-a + √(-24 + a^2))
 (1//2)*(-a - √(-24 + a^2))
```
```jldoctest
julia> symbolic_solve(x^7 - 1, x)
2-element Vector{Any}:
  roots_of((1//1) + x + x^2 + x^3 + x^4 + x^5 + x^6, x)
 1
```
### `solve_multivar` (uses Groebner basis and `solve_univar` to find roots)
!!! note
    Similar to `solve_univar`, `Groebner` is needed for `solve_multivar` or to be fully functional.
```jldoctest
julia> using Groebner

julia> @variables x y z
3-element Vector{Num}:
 x
 y
 z

julia> eqs = [x+y^2+z, z*x*y, z+3x+y]
3-element Vector{Num}:
 x + z + y^2
       x*y*z
  3x + y + z

julia> symbolic_solve(eqs, [x,y,z])
3-element Vector{Any}:
 Dict{Num, Any}(z => 0, y => 1//3, x => -1//9)
 Dict{Num, Any}(z => 0, y => 0, x => 0)
 Dict{Num, Any}(z => -1, y => 1, x => 0)
```

!!! note
    If `Nemo` or `Groebner` are not imported when needed, the solver throws an error.
```jldoctest
julia> using Symbolics

julia> @variables x y z;

julia> symbolic_solve(x+1, x)
ERROR: "Nemo is required. Execute `using Nemo` to enable this functionality."

julia> symbolic_solve([x+1, y], [x, y])
ERROR: "Groebner bases engine is required. Execute `using Groebner` to enable this functionality."
```
### `solve_multipoly` (uses GCD between the input polys)
```jldoctest
julia> symbolic_solve([x-1, x^3 - 1, x^2 - 1, (x-1)^20], x)
1-element Vector{BigInt}:
 1
```

### `ia_solve` (solving by isolation and attraction)
```jldoctest
julia> symbolic_solve(2^(x+1) + 5^(x+3), x)
1-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
 (-slog(2) - log(complex(-1)) + 3slog(5)) / (slog(2) - slog(5))
```
```jldoctest
julia> symbolic_solve(log(x+1)+log(x-1), x)
2-element Vector{SymbolicUtils.BasicSymbolic{BigFloat}}:
 (1//2)*√(8.0)
 (-1//2)*√(8.0)
```
```jldoctest
julia> symbolic_solve(a*x^b + c, x)
((-c)^(1 / b)) / (a^(1 / b))
```

### Evaluating output (converting to floats)
If you want to evaluate the exact expressions found by `symbolic_solve`, you can do the following:
```jldoctest
julia> roots = symbolic_solve(2^(x+1) + 5^(x+3), x)
1-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
 (-slog(2) - log(complex(-1)) + 3slog(5)) / (slog(2) - slog(5))

julia> eval.(Symbolics.toexpr.(roots))
1-element Vector{Complex{BigFloat}}:
 -4.512941594732059759689023145584186058252768936052415430071569066192919491762214 + 3.428598090438030380369414618548038962770087500755160535832807433942464545729382im
```
"""
function symbolic_solve(expr, x::T; dropmultiplicity = true, warns = true) where {T}
    expr_univar = false
    x_univar = false

    if (T == Num || T == SymbolicUtils.BasicSymbolic{Real})
        x_univar = true
        check_x(x)
    else
        for var in x
            check_x(var)
        end
        if length(x) == 1
            x = x[1]
            x_univar = true
        end
    end

    if !(expr isa Vector)
        expr_univar = true
        expr = expr isa Equation ? expr.lhs - expr.rhs : expr
        check_expr_validity(expr)
        isequal(expr, 0) && return []
    else
        expr = Vector{Any}(expr)
        for i in eachindex(expr)
            expr[i] = expr[i] isa Equation ? expr[i].lhs - expr[i].rhs : expr[i]
            check_expr_validity(expr[i])
            if !check_poly_inunivar(expr[i], x)
                warns && @warn("Solve can not solve this input currently")
                return nothing
            end
        end
        expr = Vector{Num}(expr)
    end

    if expr_univar && !x_univar
        expr = [expr]
        expr_univar = false
    end

    if x_univar
        sols = []
        if expr_univar
            sols = check_poly_inunivar(expr, x) ?
                   solve_univar(expr, x, dropmultiplicity = dropmultiplicity) :
                   ia_solve(expr, x, warns = warns)
            isequal(sols, nothing) && return nothing
        else
            for i in eachindex(expr)
                if !check_poly_inunivar(expr[i], x)
                    warns && @warn("Solve can not solve this input currently")
                    return nothing
                end
            end
            sols = solve_multipoly(
                expr, x, dropmultiplicity = dropmultiplicity, warns = warns)
            isequal(sols, nothing) && return nothing
        end

        sols = map(postprocess_root, sols)
        return sols
    end

    if !expr_univar && !x_univar
        for e in expr
            for var in x
                if !check_poly_inunivar(e, var)
                    warns && @warn("This system can not be currently solved by solve.")
                    return nothing
                end
            end
        end

        sols = solve_multivar(expr, x, dropmultiplicity = dropmultiplicity)
        isequal(sols, nothing) && return nothing
        for sol in sols
            for var in x
                sol[var] = postprocess_root(sol[var])
            end
        end

        return sols
    end
end

function symbolic_solve(expr; x...)
    r = filter_poly.(expr)
    subs, filtered = r isa Tuple ? r : (map(t -> t[1], r), map(t -> t[2], r))

    sub_vars = []
    if !isempty(subs)
        sub_vars = reduce(vcat, subs isa Dict ? collect(keys(subs)) : collect.(keys.(subs)))
    end

    vars_list = get_variables.(filtered)
    filt_vars = unique(isa(vars_list, AbstractArray) && all(isa.(vars_list, AbstractArray)) ? reduce(vcat, vars_list) : vars_list)

    vars = isempty(sub_vars) ? filt_vars : setdiff(filt_vars, sub_vars)
    vars = wrap.(vars)
    @assert all(v isa Num for v in vars) "All variables should be Nums or BasicSymbolics"

    return symbolic_solve(expr, vars; x...)
end

"""
    solve_univar(expression, x; dropmultiplicity=true)
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

- dropmultiplicity (optional): Print repeated roots or not?

# Examples

"""
function solve_univar(expression, x; dropmultiplicity = true)
    args = []
    mult_n = 1
    expression = unwrap(expression)
    expression = expression isa PolyForm ? SymbolicUtils.toterm(expression) : expression

    # handle multiplicities (repeated roots), i.e. (x+1)^20
    if iscall(expression)
        expr = unwrap(simplify((copy(wrap(expression)))))
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

        # multiplicities (repeated roots)
        if !dropmultiplicity
            og_arr_roots = copy(arr_roots)
            for i in 1:(mult_n - 1)
                append!(arr_roots, og_arr_roots)
            end
        end

        return arr_roots
    end

    if length(factors) != 1
        for factor in factors_subbed
            roots = solve_univar(factor, x, dropmultiplicity = dropmultiplicity)
            append!(arr_roots, roots)
        end
    end

    if isequal(arr_roots, [])
        return [RootsOf(wrap(expression), wrap(x))]
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
function solve_multipoly(polys::Vector, x::Num; dropmultiplicity = true, warns = true)
    polys = unique(polys)

    if length(polys) < 1
        warns && @warn("No expressions entered")
        return nothing
    end
    if length(polys) == 1
        return solve_univar(polys[1], x, dropmultiplicity = dropmultiplicity)
    end

    gcd = gcd_use_nemo(polys[1], polys[2])

    for i in eachindex(polys)[3:end]
        gcd = gcd_use_nemo(gcd, polys[i])
    end

    if isequal(gcd, 1)
        return []
    end

    return solve_univar(gcd, x, dropmultiplicity = dropmultiplicity)
end

function solve_multivar(eqs::Any, vars::Any; dropmultiplicity = true, warns = true)
    throw("Groebner bases engine is required. Execute `using Groebner` to enable this functionality.")
end
