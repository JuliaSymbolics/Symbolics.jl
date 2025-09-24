module SymbolicsGroebnerExt

using Groebner
const Nemo = Groebner.Nemo
using Symbolics
using Symbolics: Num, symtype, BasicSymbolic
import Symbolics.PrecompileTools
import Symbolics.Bijections
import Symbolics: SymbolicUtils
import Symbolics: DP, MP

function Symbolics.groebner_basis(polynomials::Vector{Num}; ordering=InputOrdering(), kwargs...)
    polynoms, poly_to_bs = Symbolics.symbol_to_poly(polynomials)
    basis = groebner_basis_poly(polynoms, poly_to_bs; ordering, kwargs...)
    Symbolics.poly_to_symbol(basis, poly_to_bs)
end

function groebner_basis_poly(polynoms::Vector{<:DP.Polynomial}, poly_to_bs::Bijections.Bijection; ordering=InputOrdering(), kwargs...)
    bs_to_poly = Bijections.active_inv(poly_to_bs)
    ordering = Groebner.ordering_transform(ordering, bs_to_poly)
    return Groebner.groebner(polynoms; ordering=ordering, kwargs...)
end

"""
    is_groebner_basis(polynomials; kwargs...)

Checks whether the given `polynomials` forms a Groebner basis using Groebner.jl
as the backend.

## Optional Arguments

The Groebner.jl backend provides a number of useful keyword arguments, which are
also available for this function. See `?Groebner.isgroebner`.

## Example

```jldoctest
julia> using Symbolics, Groebner

julia> @variables x y;

julia> is_groebner_basis([x^2 - y^2, x*y^2 + x, y^3 + y])
```
"""
function Symbolics.is_groebner_basis(polynomials::Vector{<:Union{Num, BasicSymbolic}}; kwargs...)
    polynoms, _ = Symbolics.symbol_to_poly(polynomials)
    Groebner.isgroebner(polynoms; kwargs...)
end

### Solver ###

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
    if num isa Num || den isa Num
        num / den
    else
        num // den
    end
end

function nemo_crude_evaluate(poly::Nemo.ZZRingElem, varmap)
    BigInt(poly)
end

function gen_separating_var(vars)
    n = 1
    new_var = (Symbolics.@variables _T)[1]
    present = any(isequal(new_var, var) for var in vars)
    while present
        new_var = Symbolics.variables(repeat("_", n) * "_T")[1]
        present = any(isequal(new_var, var) for var in vars)
        n += 1
    end
    return new_var
end

# Given a GB in k[params][vars] produces a GB in k(params)[vars]
function demote(gb, gb_as_poly, bs_to_poly, vars::Vector{Num}, params::Vector{Num})
    length(gb) == 1 && SymbolicUtils._isone(gb[1]) && return gb
    isequal(gb, [1]) && return gb
    # gb = Symbolics.wrap.(SymbolicUtils.toterm.(gb))
    Symbolics.check_polynomial.(gb)

    all_vars = [vars..., params...]
    nemo_ring, nemo_all_vars = Nemo.polynomial_ring(Nemo.QQ, map(string, all_vars))

    sym_to_nemo = Dict(all_vars .=> nemo_all_vars)
    nemo_to_sym = Dict(v => k for (k, v) in sym_to_nemo)
    pvar_to_nemo = [bs_to_poly[v] for v in all_vars] => nemo_all_vars
    nemo_gb = [poly(pvar_to_nemo) for poly in gb_as_poly]

    nemo_vars = view(nemo_all_vars, 1:length(vars))
    nemo_params = view(nemo_all_vars, length(vars)+1:length(all_vars))

    ring_flat = parent(nemo_vars[1])
    ring_param, params_demoted = Nemo.polynomial_ring(Nemo.base_ring(ring_flat), map(string, nemo_params))
    ring_demoted, vars_demoted = Nemo.polynomial_ring(Nemo.fraction_field(ring_param), map(string, nemo_vars), internal_ordering=:lex)
    varmap = Dict(vcat(nemo_vars .=> vars_demoted, nemo_params .=> params_demoted))
    gb_demoted = map(f -> ring_demoted(nemo_crude_evaluate(f, varmap)), nemo_gb)
    result = empty(gb_demoted)
    while true
        gb_demoted = map(f -> Nemo.map_coefficients(c -> c // Nemo.leading_coefficient(f), f), gb_demoted)
        for i in 1:length(gb_demoted)
            f = gb_demoted[i]
            f_nf = Nemo.normal_form(f, result)
            if !iszero(f_nf)
                push!(result, f_nf)
            end
        end
        isequal(gb_demoted, result) && break
        gb_demoted = result
        result = empty(result)
    end
    @assert all(f -> isone(Nemo.leading_coefficient(f)), result)

    sym_to_nemo = Dict(all_vars .=> [vars_demoted; params_demoted])
    nemo_to_sym = Dict(v => k for (k, v) in sym_to_nemo)

    final_result = Num[]

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

function solve_zerodim(eqs::Vector, vars::Vector{Num}; dropmultiplicity=true, warns=true)
    # Reference: Rouillier, F. Solving Zero-Dimensional Systems
    # Through the Rational Univariate Representation.
    # AAECC 9, 433–461 (1999). https://doi.org/10.1007/s002000050114
    
    rng = Groebner.Random.Xoshiro(42)
    all_indeterminates = collect(reduce(union!, map(Symbolics.get_variables, eqs)))
    params = map(Symbolics.Num ∘ Symbolics.wrap, setdiff(all_indeterminates, vars))

    # Use a new variable to separate the input polynomials (Reference above)
    new_var = gen_separating_var(vars)
    old_len = length(vars)
    old_vars = deepcopy(vars)
    vars = vcat(vars, new_var)

    new_eqs = []
    generating = true
    n_iterations = 1
    separating_form = new_var
    eqs = Symbolics.wrap.(eqs)

    while generating
        new_eqs = copy(eqs)
        separating_form = new_var
        for i = 1:(old_len)
            separating_form += BigInt(rand(rng, -n_iterations:n_iterations))*vars[i]
        end

        if isequal(separating_form, new_var)
            continue
        end

        push!(new_eqs, separating_form)

        polynoms, poly_to_bs = Symbolics.symbol_to_poly(new_eqs)
        basis = groebner_basis_poly(polynoms, poly_to_bs; ordering = Lex(vcat(vars, params)))
        new_eqs = Symbolics.poly_to_symbol(basis, poly_to_bs)

        # handle "unsolvable" case
        if SymbolicUtils._iszero(new_eqs[1])
            return []
        end

        for i in reverse(eachindex(new_eqs))
            all_present = Symbolics.get_variables(new_eqs[i])
            if length(intersect(all_present, vars)) < 1
                deleteat!(new_eqs, i)
                deleteat!(basis, i)
            end
        end
        new_eqs = demote(new_eqs, basis, Bijections.active_inv(poly_to_bs), vars, params)
        new_eqs = map(Symbolics.unwrap, new_eqs)
        # condition for positive dimensionality, i.e. infinite solutions
        if length(new_eqs) < length(vars)
            warns && @warn("Infinite number of solutions")
            return nothing
        end

        # Exit in the Shape Lemma case:
        # g(T, params) = 0
        # x1 - f1(T, params) = 0
        # ...
        # xn - fn(T, params) = 0
        generating = !(length(new_eqs) == length(vars))
        if length(new_eqs) == length(vars)
            vars_in_1 = Symbolics.get_variables(new_eqs[1])
            setdiff!(vars_in_1, params)
            generating |= !(length(vars_in_1) == 1 && isequal(first(vars_in_1), new_var))
            for i in eachindex(new_eqs)[2:end]
                present_vars = setdiff(Symbolics.get_variables(new_eqs[i]), new_var)
                present_vars = setdiff(present_vars, params)
                isempty(present_vars) && (generating = false; break;)
                var_i = first(present_vars)
                condition1 = length(present_vars) == 1
                condition2 = Symbolics.degree(new_eqs[i], var_i) == 1
                generating |= !(condition1 && condition2)
            end
        end

        # non-cyclic case
        if n_iterations > 10 
            warns && @warn("symbolic_solve can not currently solve this system of polynomials.")
            return nothing
        end

        n_iterations += 1
    end

    solutions = []

    # first, solve the first minimal polynomial
    @assert length(new_eqs) == length(vars)
    vars_in_1 = Symbolics.get_variables(new_eqs[1])
    setdiff!(vars_in_1, params)
    @assert length(vars_in_1) == 1 && isequal(first(vars_in_1), new_var)
    minpoly_sols = Symbolics.symbolic_solve(Symbolics.wrap(new_eqs[1]), new_var, dropmultiplicity=dropmultiplicity)
    solutions = [Dict{Num, Any}(new_var => sol) for sol in minpoly_sols]

    new_eqs = new_eqs[2:end]

    # second, iterate over eqs and sub each found solution
    # then add the roots of the remaining unknown variables 
    for (i, eq) in enumerate(new_eqs)
        present_vars = setdiff!(Symbolics.get_variables(eq), params)
        present_vars = setdiff!(present_vars, new_var)
        @assert length(present_vars) == 1
        var_tosolve = first(present_vars)
        @assert Symbolics.degree(eq, var_tosolve) == 1
        @assert !isempty(solutions)
        for roots in solutions
            subbded_eq = Symbolics.substitute(eq, Dict([new_var => roots[new_var]]); fold=false)
            subbded_eq = Symbolics.substitute(subbded_eq, Dict([var_tosolve => 0]); fold=false)
            new_var_sols = [-subbded_eq]
            @assert length(new_var_sols) == 1
            root = new_var_sols[1]
            roots[var_tosolve] = root
        end
    end

    vars = vars[1:end-1]
    for roots in solutions
        delete!(roots, new_var)
    end

    return solutions
end

function transendence_basis(sys, vars)
    J = Symbolics.jacobian(sys, vars)
    x0 = Dict(v => rand(-10:10) for v in vars)
    J_x0 = map(Symbolics.value, substitute(J, x0))
    rk, rref = Nemo.rref(Nemo.matrix(Nemo.QQ, J_x0))
    pivots = Int[]
    for i in 1:length(sys)
        col = findfirst(!iszero, rref[i, :])
        !isnothing(col) && push!(pivots, col)
    end
    vars[setdiff(collect(1:length(vars)), pivots)]
end

function Symbolics.solve_multivar(eqs::Vector, vars::Vector{Num}; dropmultiplicity=true, warns=true)
    sol = solve_zerodim(eqs, vars; dropmultiplicity=dropmultiplicity, warns=warns)
    !isnothing(sol) && return sol
    tr_basis = transendence_basis(eqs, vars)
    isempty(tr_basis) && return nothing
    vars_gen = setdiff(vars, tr_basis)
    sol = solve_zerodim(eqs, vars_gen; dropmultiplicity=dropmultiplicity, warns=warns)
    sol === nothing && return nothing
    for roots in sol
        for x in tr_basis
            roots[x] = x
        end
    end

    sol
end

end # module
