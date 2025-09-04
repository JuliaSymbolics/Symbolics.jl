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

    return _get_der_order(fast_substitute(expr, Dict(Differential(t)(x) => x)), x, t) + 1
end

function reduce_rule(expr, Dt)
    iscall(expr) && isequal(operation(expr), Dt) ? wrap(arguments(expr)[1]) : nothing
end

"""
    unwrap_der(expr, Dt)

Helper function to unwrap derivatives of `f(t)` in `expr` with respect to the differential operator `Dt = Differential(t)`. Returns a tuple `(n, base_expr)`, where `n` is the order of the derivative and `base_expr` is the expression with the derivatives removed. If `expr` does not contain `f(t)` or its derivatives, returns `(0, expr)`.
"""
function unwrap_der(expr, Dt)

    if reduce_rule(unwrap(expr), Dt) === nothing
        return 0, expr
    end

    order, expr = unwrap_der(reduce_rule(unwrap(expr), Dt), Dt)
    return order + 1, expr
end

# takes into account fractions
function _true_factors(expr)
    expr = flatten_fractions(unwrap(expr)) # flatten nested fractions

    numerator_factors = SymbolicUtils.numerators(unwrap(expr))
    denominator_factors = SymbolicUtils.denominators(unwrap(expr))

    facs = filter(fac -> !isequal(fac, 1), [numerator_factors; 1 ./ denominator_factors])
    return isempty(facs) ? [1] : facs
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
    eq = fast_substitute(eq, y_sub)
    
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

    return fast_substitute(expr, rev_y_sub)
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
    if iscall(expr) && isequal(operation(expr), sin) && any(isequal.(t, factors(arguments(expr)[1])))
        return arguments(expr)[1]/t, true
    end

    if iscall(expr) && isequal(operation(expr), cos) && any(isequal.(t, factors(arguments(expr)[1])))
        return arguments(expr)[1]/t, false
    end

    return nothing
end