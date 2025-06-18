# Checks that the expression is a polynomial with integer or rational
# coefficients
function check_polynomial(poly; strict=true)
    poly = wrap(poly)
    vars = get_variables(poly)
    distr, rem = polynomial_coeffs(poly, vars)
    if strict
        @assert isequal(rem, 0) "Not a polynomial"
        @assert all(c -> c isa Integer || c isa Rational, collect(values(distr))) "Coefficients must be integer or rational"
        return true
    else
        return isequal(rem, 0)
    end
end

# factor(x^2*y + b*x*y - a*x - a*b)  ->  (x*y - a)*(x + b)
function factor_use_nemo(poly::Any)
    throw("Nemo is required. Execute `using Nemo` to enable this functionality.")
end

