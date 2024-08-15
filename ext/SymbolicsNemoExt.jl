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


function Symbolics.demote(gb, vars::Vector{Num}, params::Vector{Num})
    gb = Symbolics.wrap.(SymbolicUtils.toterm.(gb))
    Symbolics.check_polynomial.(gb)

    all_vars = [vars..., params...]
    nemo_ring, nemo_all_vars = Nemo.polynomial_ring(Nemo.QQ, map(string, all_vars))

    sym_to_nemo = Dict(all_vars .=> nemo_all_vars)
    nemo_to_sym = Dict(v => k for (k, v) in sym_to_nemo)
    nemo_gb = Symbolics.substitute(gb, sym_to_nemo)
    nemo_gb = Symbolics.substitute(nemo_gb, sym_to_nemo)

    nemo_vars = [v for (k, v) in sym_to_nemo if any(isequal(k, var) for var in vars)]
    nemo_params = [v for (k, v) in sym_to_nemo if any(isequal(k, param) for param in params)]

    ring_flat = parent(nemo_vars[1])
    ring_param, params_demoted = Nemo.polynomial_ring(base_ring(ring_flat), map(string, nemo_params))
    ring_demoted, vars_demoted = Nemo.polynomial_ring(fraction_field(ring_param), map(string, nemo_vars), internal_ordering=Nemo.internal_ordering(ring_flat))
    varmap = Dict((nemo_vars .=> vars_demoted)..., (nemo_params .=> params_demoted)...)
    gb_demoted = map(f -> nemo_crude_evaluate(f, varmap), nemo_gb)
    result = empty(gb_demoted)
    for i in 1:length(gb_demoted)
        gb_demoted = map(f -> map_coefficients(c -> c // leading_coefficient(f), f), gb_demoted)
        f = gb_demoted[i]
        f_nf = Nemo.normal_form(f, result)
        if !iszero(f_nf)
            push!(result, f_nf)
        end
    end

    sym_to_nemo = Dict(sym => nem for sym in all_vars for nem in [vars_demoted..., params_demoted...] if isequal(string(sym),string(nem)))
    nemo_to_sym = Dict(v => k for (k, v) in sym_to_nemo)

    final_result = []

    for i in eachindex(result)

        monoms = collect(Nemo.monomials(result[i]))
        coeffs = collect(Nemo.coefficients(result[i]))

        poly = 0
        for j in eachindex(monoms)
            poly += nemo_crude_evaluate(coeffs[j], nemo_to_sym) * nemo_crude_evaluate(monoms[j], nemo_to_sym)
        end
        push!(final_result, poly)
    end
        
    final_result
end

end # module
