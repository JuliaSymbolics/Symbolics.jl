# Checks that the expression is a polynomial with integer or rational
# coefficients
function check_polynomial(poly)
    poly = wrap(poly)
    vars = get_variables(poly)
    distr, rem = polynomial_coeffs(poly, vars)
    @assert isequal(rem, 0) "Not a polynomial"
    @assert all(c -> c isa Integer || c isa Rational, collect(values(distr))) "Coefficients must be integer or rational"
    return true
end

# factor(x^2*y + b*x*y - a*x - a*b)  ->  (x*y - a)*(x + b)
function factor_use_nemo(poly::Any)
    throw("Nemo is required. Execute `using Nemo` to enable this functionality.")
end

# gcd(x^2 - y^2, x^3 - y^3) -> x - y
function gcd_use_nemo(poly1::Any, poly2::Any)
    throw("Nemo is required. Execute `using Nemo` to enable this functionality.")
end
