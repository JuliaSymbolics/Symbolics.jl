# We need this
import Nemo


# Map each variable of the given poly.
# Can be used to transform Nemo polynomial to expression.
function nemo_crude_evaluate(poly::Nemo.MPolyRingElem, varmap)
    @assert Nemo.coefficient_ring(poly) in (Nemo.ZZ, Nemo.QQ)
    new_poly = 0
    for (i, term) in enumerate(Nemo.terms(poly))
        new_term = Rational(Nemo.coeff(poly, i))
        for var in Nemo.vars(term)
            exp = Nemo.degree(term, var)
            exp == 0 && continue
            new_var = varmap[var]
            new_term *= new_var^exp
        end
        new_poly += new_term
    end
    new_poly
end

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
function factor_use_nemo(poly::Num)
    check_polynomial(poly)
    degree(poly) == 0 && return poly, Num[]
    vars = get_variables(poly)
    nemo_ring, nemo_vars = Nemo.polynomial_ring(Nemo.QQ, map(string, vars))
    sym_to_nemo = Dict(vars .=> nemo_vars)
    nemo_to_sym = Dict(v => k for (k, v) in sym_to_nemo)
    nemo_poly = substitute(poly, sym_to_nemo)
    nemo_fac = Nemo.factor(nemo_poly)
    nemo_unit = Nemo.unit(nemo_fac)
    nemo_factors = collect(keys(nemo_fac.fac)) 
    sym_unit = Rational(Nemo.coeff(nemo_unit, 1))
    sym_factors = map(f -> wrap(nemo_crude_evaluate(f, nemo_to_sym)), nemo_factors)

    for (i, fac) in enumerate(sym_factors)
        sym_factors[i] = fac^(collect(values(nemo_fac.fac))[i])
    end

    return sym_unit, sym_factors
end

find_nemo_var(symbolics_var, nemo_vars) = nemo_vars[findfirst(s -> string(s) == string(symbolics_var), nemo_vars)]

function split_by_variable(f, var)
    @assert var in Nemo.gens(Nemo.parent(f))
    F, G = zero(f), []
    for t in Nemo.terms(f)
        d = Nemo.degree(t, var)
        if d > 0
            push!(G, (d, Nemo.divexact(t, var^d)))
        else
            F += t
        end
    end
    return F,G
end


# gcd(x^2 - y^2, x^3 - y^3) -> x - y
function gcd_use_nemo(poly1::Num, poly2::Num)
    check_polynomial(poly1)
    check_polynomial(poly2)
    vars1 = get_variables(poly1)
    vars2 = get_variables(poly2)
    vars = vcat(vars1, vars2)
    nemo_ring, nemo_vars = Nemo.polynomial_ring(Nemo.QQ, map(string, vars))
    sym_to_nemo = Dict(vars .=> nemo_vars)
    nemo_to_sym = Dict(v => k for (k, v) in sym_to_nemo)
    nemo_poly1 = substitute(poly1, sym_to_nemo)
    nemo_poly2 = substitute(poly2, sym_to_nemo)
    nemo_gcd = Nemo.gcd(nemo_poly1, nemo_poly2)
    sym_gcd = wrap(nemo_crude_evaluate(nemo_gcd, nemo_to_sym))
    return sym_gcd
end
