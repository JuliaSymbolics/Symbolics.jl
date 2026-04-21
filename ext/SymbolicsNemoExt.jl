module SymbolicsNemoExt
using Nemo
import Symbolics.PrecompileTools
import Symbolics.Bijections

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
    num = nemo_crude_evaluate(numerator(poly), varmap)
    den = nemo_crude_evaluate(denominator(poly), varmap)
    # Collapse integer-valued fractions (n//1 → n) so returned expressions
    # don't carry noisy `1//1` constants — purely a presentation concern;
    # mathematically identical.
    if den isa Integer && den == 1
        return num
    end
    num // den
end

# Nemo's QQFieldElem is the concrete rational type returned as coefficients
# of polynomials over ℚ. Without this specific method, the FracElem dispatch
# above would still handle it, but going through `QQFieldElem` directly is
# cheaper and makes the narrowing explicit.
function nemo_crude_evaluate(poly::Nemo.QQFieldElem, varmap)
    n = BigInt(numerator(poly))
    d = BigInt(denominator(poly))
    if d == 1
        return typemin(Int) <= n <= typemax(Int) ? Int(n) : n
    end
    return Rational{BigInt}(n, d)
end

function nemo_crude_evaluate(poly::Nemo.ZZRingElem, varmap)
    # Nemo's integer elements were historically wrapped in a Rational; that
    # leaked into downstream expressions as `n//1` constants which display
    # awkwardly. Return a plain Int when it fits, BigInt when it doesn't —
    # both are mathematically identical to the Rational form but render
    # naturally.
    n = BigInt(poly)
    return typemin(Int) <= n <= typemax(Int) ? Int(n) : n
end

# factor(x^2*y + b*x*y - a*x - a*b)  ->  (x*y - a)*(x + b)
function Symbolics.factor_use_nemo(poly::Num)
    Symbolics.check_polynomial(poly)
    Symbolics.degree(poly) == 0 && return poly, Num[]
    mp_polys, poly_to_bs = Symbolics.symbol_to_poly([poly])
    mp_poly = only(mp_polys)
    vars = collect(Symbolics.get_variables(poly))
    bs_to_poly = Bijections.active_inv(poly_to_bs)
    poly_vars = map(Base.Fix1(getindex, bs_to_poly), vars)
    nemo_ring, nemo_vars = Nemo.polynomial_ring(Nemo.QQ, map(string, vars))
    nemo_poly = mp_poly(poly_vars => nemo_vars)
    nemo_fac = Nemo.factor(nemo_poly)
    nemo_unit = Nemo.unit(nemo_fac)
    nemo_factors = collect(keys(nemo_fac.fac)) 
    sym_unit = Rational(Nemo.coeff(nemo_unit, 1))
    nemo_to_sym = Dict(nemo_vars .=> vars)
    sym_factors = map(f -> Symbolics.wrap(nemo_crude_evaluate(f, nemo_to_sym)), nemo_factors)
    for (i, fac) in enumerate(sym_factors)
        sym_factors[i] = fac^(collect(values(nemo_fac.fac))[i])
    end

    return sym_unit, sym_factors
end

# Nemo-backed implementation of the Symbolics.factor extension hook. Called by
# `Symbolics.factor` when the pure-Julia core leaves a residual it cannot
# prove irreducible (typically degree ≥ 4 with no rational root). Returns a
# `Vector{Tuple{Num, Int}}` of (factor, multiplicity) pairs whose product
# equals `residual` — or `nothing` if Nemo fails or the input isn't monic.
function Symbolics._factor_with_nemo(residual::Num, var)
    try
        Symbolics.check_polynomial(residual)
        mp_polys, poly_to_bs = Symbolics.symbol_to_poly([residual])
        mp_poly = only(mp_polys)
        vars = collect(Symbolics.get_variables(residual))
        isempty(vars) && return nothing
        bs_to_poly = Bijections.active_inv(poly_to_bs)
        poly_vars = map(Base.Fix1(getindex, bs_to_poly), vars)
        _, nemo_vars = Nemo.polynomial_ring(Nemo.QQ, map(string, vars))
        nemo_poly = mp_poly(poly_vars => nemo_vars)
        nemo_fac = Nemo.factor(nemo_poly)

        # For a monic residual (what Symbolics.factor feeds us), Nemo's unit
        # should be the constant polynomial 1. Bail on anything else so the
        # caller falls back to returning the residual unsplit.
        nemo_unit = Nemo.unit(nemo_fac)
        unit_const = Rational(Nemo.coeff(nemo_unit, 1))
        unit_const == 1 || return nothing

        nemo_to_sym = Dict(nemo_vars .=> vars)
        result = Tuple{Num, Int}[]
        for (fac, mult) in nemo_fac.fac
            p = Symbolics.wrap(nemo_crude_evaluate(fac, nemo_to_sym))
            push!(result, (p, Int(mult)))
        end
        return result
    catch
        return nothing
    end
end

# Helps with precompilation time
# PrecompileTools.@setup_workload begin
#     @variables a b c x y z
#     expr_with_params = expand((x + b)*(x^2 + 2x + 1)*(x^2 - a))
#     equation1 = a*log(x)^b + c ~ 0
#     equation_polynomial = 9^x + 3^x + 2
#     exp_eq = 5*2^(x+1) + 7^(x+3)
#     PrecompileTools.@compile_workload begin
#         symbolic_solve(equation1, x)
#         symbolic_solve(equation_polynomial, x)
#         symbolic_solve(exp_eq)
#         symbolic_solve(expr_with_params, x, dropmultiplicity=false)
#         symbolic_solve(x^10 - a^10, x, dropmultiplicity=false)
#     end
# end

end # module
