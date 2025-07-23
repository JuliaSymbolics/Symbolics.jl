import DomainSets.ClosedInterval

function laplace(expr, f, t, s, F)
    Dt = Differential(t)
    Ds = Differential(s)
    # from https://tutorial.math.lamar.edu/Classes/DE/Laplace_Table.aspx
    transform_rules = Symbolics.Chain([
        @rule 1 => 1/s
        @rule exp(t) => 1/(s - 1)
        @rule exp(~a * t) => 1/(-~a + s)
        @rule t => 1/s^2
        @rule t^~n => factorial(~n)/s^(~n + 1)
        @rule sqrt(t) => sqrt(pi)/(2 * s^(3/2))
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
        @rule ~g * exp(~c*t) => laplace(~g, f, t, s - ~c, F) # s-shift rule
        @rule t*f(t) => -Ds(F(s)) # s-derivative rule
        @rule t^(~n)*f(t) => (-1)^(~n) * (Ds^~n)(F(s)) # s-derivative rule
        @rule f(~a + t) => exp(~a*s)*F(s) # t-shift rule
        @rule f(t) => F(s)
    ])

    transformed = transform_rules(expr)
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
        factors = _true_factors(term)
        constant = filter(x -> isempty(Symbolics.get_variables(x)), factors)
        if !isempty(constant)
            result += laplace(term / constant[1], f, t, s, F) * constant[1]
        else 
            result += laplace(term, f, t, s, F)
        end
    end

    return result
end

function laplace(expr::Equation, f, t, s, F)
    return laplace(expr.lhs, f, t, s, F) ~ laplace(expr.rhs, f, t, s, F)
end

function inverse_laplace(expr, F, t, s, f)
    inverse_transform_rules = Symbolics.Chain([
        @rule 1/s => 1
        @rule 1/(~a + s) => exp(~a * t)
        @rule 1/s^(~n) => t^(~n-1) / factorial(~n-1)
        @rule 1/(2 * s^(3/2)) => sqrt(t)/sqrt(pi)
        @rule 1/((~a)^2 + s^2) => sin(~a * t)/~a
        @rule s/((~a)^2 + s^2) => cos(~a * t)
        @rule s / ((~a)^2 + s^2)^2 => t*sin(~a * t)/(2*~a)
        @rule (-(~a)^2 + s^2) / ((~a)^2 + s^2)^2 => t*cos(~a * t)
        @rule 1 / ((~a)^2 + s^2)^2 => (sin(~a*t) - ~a*t*cos(~a*t))/ (2*(~a)^3)
        @rule s^2 / ((~a)^2 + s^2)^2 => (sin(~a*t) + ~a*t*cos(~a*t)) / (2*~a)
        @rule s*((~a)^2 + s^2) / ((~a)^2 + s^2)^2 => cos(~a*t) - ~a*t*sin(~a*t)
        @rule s*(3*(~a)^2 + s^2) / ((~a)^2 + s^2)^2 => cos(~a*t) + ~a*t*sin(~a*t)
        @rule (s*sin(~b) + ~a*cos(~b)) / ((~a)^2 + s^2) => sin(~b + ~a*t)
        @rule (s*cos(~b) - ~a*sin(~b)) / ((~a)^2 + s^2) => cos(~b + ~a*t)
        @rule 1/(s^2 - (~b)^2) => sinh(~b * t)/~b
        @rule s/(s^2 - (~b)^2) => cosh(~b * t)
        @rule 1 / ((~c+s)^2 + (~b)^2) => exp(-~c*t) * sin(~b * t) / ~b
        @rule (~c+s) / ((~c+s)^2 + (~b)^2) => exp(-~c*t) * cos(~b * t)
        @rule 1 / ((~c+s)^2 - (~b)^2) => exp(-~c*t) * sinh(~b * t) / ~b
        @rule (~c+s) / ((~c+s)^2 - (~b)^2) => exp(-~c*t) * cosh(~b * t)
        @rule 1 / (~a + s)^(~n) => t^(~n-1) * exp(-~a * t) / factorial(~n-1)
    ])

    transformed = inverse_transform_rules(expr)
    if !isequal(transformed, expr)
        return transformed
    end

    terms = Symbolics.terms(expr)
    result = 0
    if length(terms) == 1 && length(_true_factors(terms[1])) == 1
        return f
    end
    for term in terms
        factors = _true_factors(term)
        constant = filter(x -> isempty(Symbolics.get_variables(x)), factors)
        if !isempty(constant)
            result += inverse_laplace(term / constant[1], F, t, s, f) * constant[1]
        else 
            result += inverse_laplace(term, F, t, s, f)
        end
    end

    return result
end

function inverse_laplace(expr::Equation, F, t, s, f)
    return inverse_laplace(expr.lhs, F, t, s, f) ~ inverse_laplace(expr.rhs, F, t, s, f)
end

function unwrap_der(expr, Dt)
    reduce_rule = @rule Dt(~x) => ~x

    if reduce_rule(expr) === nothing
        return 0, expr
    end

    order, expr = unwrap_der(reduce_rule(expr), Dt)
    return order + 1, expr
end