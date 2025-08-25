import DomainSets.ClosedInterval

# from https://tutorial.math.lamar.edu/Classes/DE/Laplace_Table.aspx
transform_rules(f, t, F, s) = Symbolics.Chain([
    @rule 1 => 1/s
    @rule exp(t) => 1/(s - 1)
    @rule exp(~a * t) => 1/(-~a + s)
    @rule t => 1/s^2
    @rule t^~n => factorial(~n)/s^(~n + 1)
    @rule sqrt(t) => term(sqrt, pi)/(2 * s^(3/2))
    @rule sin(t) => 1/(1 + s^2)
    @rule sin(~a * t) => ~a/((~a)^2 + s^2)
    @rule cos(t) => s/(1 + s^2)
    @rule cos(~a * t) => s/((~a)^2 + s^2)
    @rule t*sin(t) => 1/(1 + s^2)^2
    @rule t*sin(~a * t) => 2*~a*s / ((~a)^2 + s^2)^2
    @rule t*cos(t) => (s^2 - 1) / (1 + s^2)^2
    @rule t*cos(~a * t) => (-(~a)^2 + s^2) / ((~a)^2 + s^2)^2
    @rule sin(t) - t*cos(t) => 2 / (1 + s^2)^2
    @rule sin(~a*t) - ~a*t*cos(~a*t) => 2*(~a)^3 / ((~a)^2 + s^2)^2
    @rule sin(t) + t*cos(t) => 2s^2 / (1 + s^2)^2
    @rule sin(~a*t) + ~a*t*cos(~a*t) => 2*~a*s^2 / ((~a)^2 + s^2)^2
    @rule cos(~a*t) - ~a*t*sin(~a*t) => s*((~a)^2 + s^2) / ((~a)^2 + s^2)^2
    @rule cos(~a*t) + ~a*t*sin(~a*t) => s*(s^2 + 3*(~a)^2) / ((~a)^2 + s^2)^2
    @rule sin(~b + ~a*t) => (s*sin(~b) + ~a*cos(~b)) / ((~a)^2 + s^2)
    @rule cos(~b + ~a*t) => (s*cos(~b) - ~a*sin(~b)) / ((~a)^2 + s^2)
    @rule sinh(~a * t) => ~a/(-(~a)^2 + s^2)
    @rule cosh(~a * t) => s/(-(~a)^2 + s^2)
    @rule exp(~a*t) * sin(~b * t) => ~b / ((~b)^2 + (-~a+s)^2)
    @rule exp(~a*t) * cos(~b * t) => (-~a+s) / ((~b)^2 + (-~a+s)^2)
    @rule exp(~a*t) * sinh(~b * t) => ~b / (-(~b)^2 + (-~a+s)^2)
    @rule exp(~a*t) * cosh(~b * t) => (-~a+s) / (-(~b)^2 + (-~a+s)^2)
    @rule t^~n * exp(~a * t) => factorial(~n) / (-~a + s)^(~n + 1)
    @rule t*exp(~a * t) => 1 / (-~a + s)^(2)
    @rule t^~n * exp(t) => factorial(~n) / (s)^(~n + 1)
    @rule t*exp(t) => 1 / (s)^(2)
    @rule exp(~c*t) * ~g => laplace(~g, f, t, F, s - ~c) # s-shift rule
    @rule t*f(t) => -Ds(F(s)) # s-derivative rule
    @rule t^(~n)*f(t) => (-1)^(~n) * (Ds^~n)(F(s)) # s-derivative rule
    @rule f(~a + t) => exp(~a*s)*F(s) # t-shift rule
    @rule f(t) => F(s)
])

"""
    laplace(expr, f, t, F, s)

Performs the Laplace transform of `expr` with respect to the variable `t`, where `f(t)` is a function in `expr` being transformed, and `F(s)` is the Laplace transform of `f(t)`. Returns the transformed expression in terms of `s`.

Note that `f(t)` and `F(s)` should be defined using `@syms`

Currently relies mostly on linearity and a rules table. When the rules table does not apply, it falls back to the integral definition of the Laplace transform.
"""
function laplace(expr, f, t, F, s)
    expr = expand(expr)
    Dt = Differential(t)
    Ds = Differential(s)

    transformed = transform_rules(f, t, F, s)(expr)
    if !isequal(transformed, expr)
        return transformed
    end

    # t-derivative rule
    n, expr = unwrap_der(expr, Dt)
    if n != 0 && isequal(expr, f(t))
        f0 = Symbolics.variables(:f0, 0:(n-1))
        transformed = s^n*F(s)
        for i = 1:n
            transformed -= s^(n-i)*f0[i]
        end

        return transformed
    end

    terms = Symbolics.terms(expr)
    result = 0
    if length(terms) == 1 && length(filter(x->isempty(Symbolics.get_variables(x)), _true_factors(terms[1]))) == 0
        return Integral(t in ClosedInterval(0, Inf))(expr*exp(-s*t))
    end
    for term in terms
        factors = _true_factors(wrap(term))
        constant = filter(x -> isempty(Symbolics.get_variables(x)), factors)
        if !isempty(constant)
            result += laplace(term / constant[1], f, t, F, s) * constant[1]
        else 
            result += laplace(term, f, t, F, s)
        end
    end

    return result
end

function laplace(expr::Equation, f, t, F, s)
    return laplace(expr.lhs, f, t, F, s) ~ laplace(expr.rhs, f, t, F, s)
end

# postprocess_root prevents automatic evaluation of sqrt to its floating point value
function processed_sqrt(x)
    return postprocess_root(term(sqrt, x))
end

# F and f aren't used here, but are here for future-proofing
inverse_transform_rules(F, s, f, t) = Symbolics.Chain([
    @rule 1/s => 1
    @rule 1/(~a + s) => exp(-~a * t)
    @rule 1/s^(~n) => t^(~n-1) / factorial(~n-1)
    @rule 1/(2 * s^(3/2)) => sqrt(t)/term(term(sqrt, pi))
    @rule 1/(~a + s^2) => sin(processed_sqrt(~a) * t)/processed_sqrt(~a)
    @rule s/(~a + s^2) => cos(processed_sqrt(~a) * t)
    @rule s / (~a + s^2)^2 => t*sin(processed_sqrt(~a) * t)/(2*processed_sqrt(~a))
    @rule (-~a + s^2) / (~a + s^2)^2 => t*cos(processed_sqrt(~a) * t)
    @rule 1 / (~a + s^2)^2 => (sin(processed_sqrt(~a)*t) - processed_sqrt(~a)*t*cos(processed_sqrt(~a)*t))/ (2*processed_sqrt(~a)^3)
    @rule s^2 / (~a + s^2)^2 => (sin(processed_sqrt(~a)*t) + processed_sqrt(~a)*t*cos(processed_sqrt(~a)*t)) / (2*processed_sqrt(~a))
    @rule s*(~a + s^2) / (~a + s^2)^2 => cos(processed_sqrt(~a)*t) - processed_sqrt(~a)*t*sin(processed_sqrt(~a)*t)
    @rule s*(3*~a + s^2) / (~a + s^2)^2 => cos(processed_sqrt(~a)*t) + processed_sqrt(~a)*t*sin(processed_sqrt(~a)*t)
    @rule (s*sin(~b) + ~a*cos(~b)) / (~a + s^2) => sin(~b + processed_sqrt(~a)*t)
    @rule (s*cos(~b) - ~a*sin(~b)) / ((~a)^2 + s^2) => cos(~b + ~a*t)
    @rule 1/(s^2 - (~b)^2) => sinh(~b * t)/~b
    @rule s/(s^2 - (~b)^2) => cosh(~b * t)
    @rule 1 / ((~c+s)^2 + (~b)^2) => exp(-~c*t) * sin(~b * t) / ~b
    @rule (~c+s) / ((~c+s)^2 + (~b)^2) => exp(-~c*t) * cos(~b * t)
    @rule 1 / ((~c+s)^2 - (~b)^2) => exp(-~c*t) * sinh(~b * t) / ~b
    @rule (~c+s) / ((~c+s)^2 - (~b)^2) => exp(-~c*t) * cosh(~b * t)
    @rule 1 / (~a + s)^(~n) => t^(~n-1) * exp(-~a * t) / factorial(~n-1)
])

"""
    inverse_laplace(expr, F, s, f, t)

Performs the inverse Laplace transform of `expr` with respect to the variable `s`, where `F(s)` is the Laplace transform of `f(t)`. Returns the transformed expression in terms of `t`.

Note that `f(t)` and `F(s)` should be defined using `@syms`.

Will perform partial fraction decomposition and linearity before applying the inverse Laplace transform rules. When unable to find a result, returns `nothing`.
"""
function inverse_laplace(expr, F, s, f, t)
    if isequal(expr, 0)
        return 0
    end
    
    # check for partial fractions
    partial_fractions = partial_frac_decomposition(expr, s)
    if partial_fractions !== nothing && !isequal(partial_fractions, expr)
        return inverse_laplace(partial_fractions, F, s, f, t)
    end

    transformed = inverse_transform_rules(F, s, f, t)(expr)
    if !isequal(transformed, expr)
        return transformed
    end

    _terms = terms(numerator(expr)) ./ denominator(expr)
    
    result = 0
    if length(_terms) == 1 && length(filter(x -> isempty(get_variables(x)), _true_factors(_terms[1]))) == 0
        @warn "Inverse laplace failed: $expr"
        return nothing # no result
    end

    # apply linearity
    for term in _terms
        factors = _true_factors(term)
        constant = filter(x -> isempty(Symbolics.get_variables(x)), factors)
        if !isempty(constant)
            result += inverse_laplace(term / constant[1], F, s, f, t) * constant[1]
        else
            result += inverse_laplace(term, F, s, f, t)
        end
    end

    return result
end

function inverse_laplace(expr::Equation, F, s, f, t)
    return inverse_laplace(expr.lhs, F, s, f, t) ~ inverse_laplace(expr.rhs, F, s, f, t)
end

"""
    laplace_solve_ode(eq, f, t, f0)
    
Solves the ordinary differential equation `eq` for the function `f(t)` using the Laplace transform method.
    
`f0` is a vector of initial conditions evaluated at `t=0` (`[f(0), f'(0), f''(0), ...]`, must be same length as order of `eq`).
"""
function laplace_solve_ode(eq, f, t, f0)
    s = variable(:𝓈)
    @syms 𝓕(s)
    transformed_eq = laplace(eq, f, t, 𝓕, s)
    transformed_eq = fast_substitute(transformed_eq, Dict(𝓕(s) => variable(:𝓕), [variable(:f0, i-1) => f0[i] for i=1:length(f0)]...))
    transformed_eq = expand(transformed_eq.lhs - transformed_eq.rhs)

    F_terms = 0
    other_terms = []
    for term in terms(transformed_eq)
        if isempty(get_variables(term, [variable(:𝓕)]))
            push!(other_terms, -1*term)
        else
            F_terms += term/variable(:𝓕) # assumes term is something times F
        end
    end

    if isempty(other_terms)
        other_terms = 0
    end

    transformed_soln = simplify(sum(other_terms ./ F_terms))

    return expand(inverse_laplace(transformed_soln, 𝓕, s, f, t))
end