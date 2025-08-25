# recursively find highest derivative order in `expr`
function _get_der_order(expr, x, t)
    if !hasderiv(unwrap(expr))
        return 0
    end

    if length(terms(expr)) > 1
        return maximum(_get_der_order.(terms(expr), Ref(x), Ref(t)))
    end

    if length(factors(expr)) > 1
        return maximum(_get_der_order.(factors(expr), Ref(x), Ref(t)))
    end

    return _get_der_order(substitute(expr, Dict(Differential(t)(x) => x)), x, t) + 1
end

"""
    unwrap_der(expr, Dt)

Helper function to unwrap derivatives of `f(t)` in `expr` with respect to the differential operator `Dt = Differential(t)`. Returns a tuple `(n, base_expr)`, where `n` is the order of the derivative and `base_expr` is the expression with the derivatives removed. If `expr` does not contain `f(t)` or its derivatives, returns `(0, expr)`.
"""
function unwrap_der(expr, Dt)
    reduce_rule = @rule Dt(~x) => ~x

    if reduce_rule(expr) === nothing
        return 0, expr
    end

    order, expr = unwrap_der(reduce_rule(expr), Dt)
    return order + 1, expr
end

# takes into account fractions
function _true_factors(expr)
    facs = factors(expr)
    true_facs::Vector{Real} = []
    frac_rule = @rule (~x)/(~y) => [~x, 1/~y]
    for fac in facs
        frac = frac_rule(fac)
        if frac !== nothing && !isequal(frac[1], 1)
            append!(true_facs, _true_factors(frac[1]))
            append!(true_facs, _true_factors(frac[2]))
        else
            push!(true_facs, fac)
        end
    end

    return convert(Vector{Num}, true_facs)
end

"""
    reduce_order(eq, x, t, ys)

Reduce order of an ODE by substituting variables for derivatives to form a system of first order ODEs
"""
function reduce_order(eq, x, t, ys)
    Dt = Differential(t)
    n = _get_der_order(eq, x, t)
    @assert n >= 1 "ODE must have at least one derivative"
    
    # reduction of order
    y_sub = Dict([[(Dt^i)(x) => ys[i+1] for i=0:n-1]; (Dt^n)(x) => variable(:𝒴)])
    eq = substitute(eq, y_sub)
    
    # isolate (Dt^n)(x)
    f = symbolic_linear_solve(eq, variable(:𝒴), check=false)
    @assert f !== nothing "Failed to isolate highest order derivative term"
    f = f[1]
    system = [ys[2:n]; f]

    return system
end

function unreduce_order(expr, x, t, ys)
    Dt = Differential(t)
    rev_y_sub = Dict(ys[i] => (Dt^(i-1))(x) for i in eachindex(ys))

    return substitute(expr, rev_y_sub)
end

function is_solution(solution, eq::Equation, x, t)
    is_solution(solution, eq.lhs - eq.rhs, x, t)
end

function is_solution(solution, eq::LinearODE)
    is_solution(solution, get_expression(eq), eq.x, eq.t)
end

function is_solution(solution, eq, x, t)
    if solution === nothing
        return false
    end

    expr = substitute(eq, Dict(x => solution))
    expr = expand(expand_derivatives(expr))
    return isequal(expr, 0)
end

function _parse_trig(expr, t)
    parse_sin = Symbolics.Chain([(@rule sin(t) => 1), (@rule sin(~x * t) => ~x)])
    parse_cos = Symbolics.Chain([(@rule cos(t) => 1), (@rule cos(~x * t) => ~x)])

    if !isequal(parse_sin(expr), expr)
        return parse_sin(expr), true
    end

    if !isequal(parse_cos(expr), expr)
        return parse_cos(expr), false
    end

    return nothing
end