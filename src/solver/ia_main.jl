using Symbolics

const SAFE_ALTERNATIVES = Dict(log => slog, sqrt => ssqrt, cbrt => scbrt)

function isolate(lhs, var; warns=true, conditions=[])
    rhs = Vector{Any}([0])
    original_lhs = deepcopy(lhs)
    lhs = unwrap(lhs)

    old_lhs = nothing
    while !isequal(lhs, var)
        subs, poly = filter_poly(lhs, var)

        if check_poly_inunivar(poly, var)
            roots = []
            new_var = gensym()
            new_var = (@variables $new_var)[1]
            lhs_roots = solve_univar(lhs - new_var, var)

            for i in eachindex(lhs_roots)
                for j in eachindex(rhs)
                    push!(roots, substitute(lhs_roots[i], Dict(new_var=>rhs[j]), fold=false))
                end
            end
            return roots, conditions
        end

        if isequal(old_lhs, lhs) 
            warns && @warn("This expression cannot be solved with the methods available to ia_solve. Try a numerical method instead.")
            return nothing
        end

        old_lhs = deepcopy(lhs)

        oper = operation(lhs)
        args = arguments(lhs)

        if oper === (+)
            for arg in args
                vars = get_variables(arg)
                if any(isequal(x, var) for x in vars)
                    continue
                end
                lhs = lhs - arg
                rhs = map(sol -> sol - arg, rhs)
            end

        elseif oper === (*)
            for arg in args
                vars = get_variables(arg)
                if any(isequal(x, var) for x in vars)
                    continue
                end
                lhs = lhs / arg
                rhs = map(sol -> sol / arg, rhs)
            end

        elseif oper === (/)
            var_innumerator = any(isequal(x, var) for x in get_variables(args[1]))
            if var_innumerator
                # x / 2 = y
                lhs = args[1]
                rhs = map(sol -> sol * args[2], rhs)
            else
                # 2 / x = y
                lhs = args[2]
                rhs = map(sol -> args[1] // sol, rhs)
            end

        elseif oper === (^)
            if any(isequal(x, var) for x in get_variables(args[1])) &&
               n_occurrences(args[2], var) == 0 && args[2] isa Integer
                lhs = args[1]
                power = args[2]
                new_roots = []

                for i in eachindex(rhs)
                    for k in 0:(args[2] - 1)
                        r = wrap(term(^, rhs[i], (1 // power)))
                        c = wrap(term(*, 2 * (k), pi)) * im / power
                        root = r * Base.MathConstants.e^c
                        push!(new_roots, root)
                    end
                end
                rhs = []
                append!(rhs, new_roots)
            elseif any(isequal(x, var) for x in get_variables(args[1])) &&
                   n_occurrences(args[2], var) == 0
                lhs = args[1]
                rhs = map(sol -> term(^, sol, 1 // args[2]), rhs)
            else
                lhs = args[2]
                rhs = map(sol -> term(/, term(slog, sol), term(slog, args[1])), rhs)
            end
        elseif has_left_inverse(oper)
            lhs = args[1]
            ia_conditions!(oper, lhs, rhs, conditions)
            invop = left_inverse(oper)
            invop = get(SAFE_ALTERNATIVES, invop, invop)
            if is_periodic(oper)
                # make this global somehow so the user doesnt need to declare it on his own
                new_var = gensym()
                new_var = (@variables $new_var)[1]
                period = fundamental_period(oper)
                rhs = map(
                    sol -> term(invop, sol) +
                           term(*, period, new_var),
                    rhs)
                @info string(new_var) * " ϵ" * " Ζ"
            else
                rhs = map(sol -> term(invop, sol), rhs)
            end
        end

        lhs = simplify(lhs)
    end

    return rhs, conditions
end

function attract(lhs, var; warns = true)
    if n_func_occ(simplify(lhs), var) <= n_func_occ(lhs, var)
        lhs = simplify(lhs)
    end
    conditions = []

    if detect_exponential(lhs, var)
        lhs = attract_exponential(lhs, var)
    end
    if detect_addlogs(lhs, var)
        lhs, new_conds = attract_logs(lhs, var)
        append!(conditions, new_conds)
    end
    lhs = attract_trig(lhs, var)

    n_func_occ(lhs, var) == 1 && return isolate(lhs, var, warns = warns, conditions=conditions)

    lhs, sub = turn_to_poly(lhs, var)

    if (isequal(sub, Dict()) || n_func_occ(lhs, collect(keys(sub))[1]) != 1)
        sqrt_poly = detect_sqrtpoly(lhs, var)
        if sqrt_poly
            return attract_and_solve_sqrtpoly(lhs, var), conditions
        else
            warns && @warn("This expression cannot be solved with the methods available to ia_solve. Try \
            a numerical method instead.")
            return nothing, conditions
        end
    end

    new_var = collect(keys(sub))[1]
    new_var_val = collect(values(sub))[1]

    roots, new_conds = isolate(lhs, new_var, warns = warns)
    append!(conditions, new_conds)
    new_roots = []

    for root in roots
        new_sol, new_conds = isolate(new_var_val - root, var, warns = warns)
        append!(conditions, new_conds)
        push!(new_roots, new_sol)
    end
    new_roots = collect(Iterators.flatten(new_roots))

    return new_roots, conditions
end

"""
    ia_solve(lhs, var)
This function attempts to solve transcendental functions by first checking
the "smart" number of occurrences in the input LHS. By smart here we mean
that polynomials are counted as 1 occurrence. for example `x^2 + 2x` is 1
occurrence of x. So we abstract all occurrences of x's as polynomials.
Say: `log(x+1) + x^2` is seen as `log(f(x)) + g(x)` so there are 2 occurrences
of x. If there is only 1 occurrence of x in an input expression, isolate is called.

Isolate reverses all operations applied on the occurrence of x until we have
`f(x) = some constant` then we can solve this using our polynomial solvers.

If more than 1 occurrence of x is found, `ia_solve` attempts to attract the
occurrences of x in order to reduce these occurrences to 1. For example,
`log(x+1) + log(x-1)` can be converted to `log(x^2 - 1)` which now could be
isolated using Isolate.

`attract(lhs, var)` currently uses 4 techniques for attraction.
- Log addition: `log(f(x)) + log(g(x)) => log(h(x))`
- Exponential simplification: `a*b^(f(x)) + c*d^(g(x)) => f(x) * log(b) - g(x) * log(d) + log(-a/c)`. And now this is actually 1 occurrence of x since `f(x)` and `g(x)` are just multiplied by constants not wrapped in some operation.
- Trig simplification: this bruteforces multiple trig identities and doesn't detect them before hand.
- Polynomialization: as a last resort, attract attempts to polynomialize the expression. Say `sin(x+2)^2 + sin(x+2) + 10` is converted to `X^2 + X + 10`, we then solve this using our polynomial solver, and afterwards, isolate `sin(x+2) = the roots found by solve for X^2 + X + 10`

After attraction, we check the number of occurrences again, and if its 1, we isolate, if not,
we throw an error to tell the user that this is currently unsolvable by our covered techniques.

# Arguments
- lhs: a Num/SymbolicUtils.BasicSymbolic
- var: variable to solve for.

# Examples
```jldoctest
julia> solve(a*x^b + c, x)
((-c)^(1 / b)) / (a^(1 / b))
```

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
julia> expr = sin(x+2)^2 + sin(x+2) + 10
10 + sin(2 + x) + sin(2 + x)^2

julia> RootFinding.ia_solve(expr, x)
[ Info: var"##230" ϵ Ζ: e.g. 0, 1, 2...
[ Info: var"##234" ϵ Ζ: e.g. 0, 1, 2...
2-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
 -2 + π*2var"##230" + asin((1//2)*(-1 + RootFinding.ssqrt(-39)))
 -2 + π*2var"##234" + asin((1//2)*(-1 - RootFinding.ssqrt(-39)))
```

All transcendental functions for which `left_inverse` is defined are supported.
To enable `ia_solve` to handle custom transcendental functions, define an inverse or
left inverse. If the function is periodic, `is_periodic` and `fundamental_period` must
be defined. If the function imposes certain conditions on its input or output (for
example, `log` requires that its input be positive) define `ia_conditions!`.

See also: [`left_inverse`](@ref), [`inverse`](@ref), [`is_periodic`](@ref),
[`fundamental_period`](@ref), [`ia_conditions!`](@ref).

# References
[^1]: [R. W. Hamming, Coding and Information Theory, ScienceDirect, 1980](https://www.sciencedirect.com/science/article/pii/S0747717189800070).
"""
function ia_solve(lhs, var; warns = true)
    nx = n_func_occ(lhs, var)
    sols = []
    conditions = []
    if nx == 0
        warns && @warn("Var not present in given expression")
        return []
    elseif nx == 1
        sols, conditions = isolate(lhs, var, warns = warns)
    elseif nx > 1
        sols, conditions = attract(lhs, var, warns = warns)
    end

    isequal(sols, nothing) && return nothing

    
    filtered_sols = []
    for i in eachindex(sols)
        if length(get_variables(sols[i])) > 0
            push!(filtered_sols, sols[i])
            continue
        end
        domain_error = false
        for j in eachindex(conditions)
            condition, t = conditions[j]
            cond_val = substitute(condition, Dict(var=>eval(toexpr(sols[i])))) 
            cond_val isa Complex && continue
            domain_error |= !t(cond_val, 0)
        end
        !domain_error && push!(filtered_sols, sols[i])
    end

    return filtered_sols
end
