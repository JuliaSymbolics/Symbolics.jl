module SymbolicsNemoExt
using Nemo
import Symbolics.PrecompileTools

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

# Helps with precompilation time
PrecompileTools.@setup_workload begin
    @variables a b c x y z
    expr_with_params = expand((x + b)*(x^2 + 2x + 1)*(x^2 - a))
    equation1 = a*log(x)^b + c ~ 0
    equation_polynomial = 9^x + 3^x + 2
    exp_eq = 5*2^(x+1) + 7^(x+3)
    PrecompileTools.@compile_workload begin
        symbolic_solve(equation1, x)
        symbolic_solve(equation_polynomial, x)
        symbolic_solve(exp_eq)
        symbolic_solve(expr_with_params, x, dropmultiplicity=false)
        symbolic_solve(x^10 - a^10, x, dropmultiplicity=false)
    end
end

end # module
