module SymbolicsNemoExt
using Nemo

if isdefined(Base, :get_extension)
    using Symbolics
    using Symbolics: Num, symtype
else
    using ..Symbolics
    using ..Symbolics: Num, symtype
end

# Map each variable of the given poly.
# Can be used to transform Nemo polynomial to expression.
function nemo_crude_evaluate(poly::Nemo.MPolyRingElem, varmap)
    new_poly = 0
    for (i, term) in enumerate(Nemo.terms(poly))
        new_term = nemo_crude_evaluate(Nemo.coeff(poly, i), varmap)
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

function nemo_crude_evaluate(poly::Nemo.FracElem, varmap)
    nemo_crude_evaluate(numerator(poly), varmap) // nemo_crude_evaluate(denominator(poly), varmap)
end

function nemo_crude_evaluate(poly::Nemo.ZZRingElem, varmap)
    Rational(poly)
end

# factor(x^2*y + b*x*y - a*x - a*b)  ->  (x*y - a)*(x + b)
function Symbolics.factor_use_nemo(poly::Num)
    Symbolics.check_polynomial(poly)
    Symbolics.degree(poly) == 0 && return poly, Num[]
    vars = Symbolics.get_variables(poly)
    nemo_ring, nemo_vars = Nemo.polynomial_ring(Nemo.QQ, map(string, vars))
    sym_to_nemo = Dict(vars .=> nemo_vars)
    nemo_to_sym = Dict(v => k for (k, v) in sym_to_nemo)
    nemo_poly = Symbolics.substitute(poly, sym_to_nemo)
    nemo_fac = Nemo.factor(nemo_poly)
    nemo_unit = Nemo.unit(nemo_fac)
    nemo_factors = collect(keys(nemo_fac.fac)) 
    sym_unit = Rational(Nemo.coeff(nemo_unit, 1))
    sym_factors = map(f -> Symbolics.wrap(nemo_crude_evaluate(f, nemo_to_sym)), nemo_factors)

    for (i, fac) in enumerate(sym_factors)
        sym_factors[i] = fac^(collect(values(nemo_fac.fac))[i])
    end

    return sym_unit, sym_factors
end

# gcd(x^2 - y^2, x^3 - y^3) -> x - y
function Symbolics.gcd_use_nemo(poly1::Num, poly2::Num)
    Symbolics.check_polynomial(poly1)
    Symbolics.check_polynomial(poly2)
    vars1 = Symbolics.get_variables(poly1)
    vars2 = Symbolics.get_variables(poly2)
    vars = vcat(vars1, vars2)
    nemo_ring, nemo_vars = Nemo.polynomial_ring(Nemo.QQ, map(string, vars))
    sym_to_nemo = Dict(vars .=> nemo_vars)
    nemo_to_sym = Dict(v => k for (k, v) in sym_to_nemo)
    nemo_poly1 = Symbolics.substitute(poly1, sym_to_nemo)
    nemo_poly2 = Symbolics.substitute(poly2, sym_to_nemo)
    nemo_gcd = Nemo.gcd(nemo_poly1, nemo_poly2)
    sym_gcd = Symbolics.wrap(nemo_crude_evaluate(nemo_gcd, nemo_to_sym))
    return sym_gcd
end

end # module
