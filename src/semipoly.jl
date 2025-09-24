import SymbolicUtils.Rewriters: RestartedChain
using DataStructures

export semipolynomial_form, semilinear_form, semiquadratic_form, polynomial_coeffs

const SemipolyDictT = Dict{BasicSymbolic{VartypeT}, BasicSymbolic{VartypeT}}

function canonicalize_poly(poly_to_bs, bs_to_poly, poly, degree, vars::AbstractSet)
    subskeys = SymbolicUtils.PolyVarT[]
    subsvals = SymbolicUtils.PolynomialT[]
    for pvar in MP.variables(poly)
        var = poly_to_bs[pvar]
        if isdiv(var)
            var = SymbolicUtils.flatten_fractions(var)
        end
        @match var begin
            BSImpl.Div(; num, den) => begin
                den_has_vars = SymbolicUtils.query(in(vars), den)
                # We only care about terms up to degree `degree`. So we only care about
                # numerators with degree <= `2degree + 1`, since if the denominator is
                # degree `degree` this entire variable goes in the residual anyway. The +1
                # is in case `degree == 0`.
                result, resid = semipolynomial_form(num, vars, 2degree + 1; consts = false)
                newpoly = SymbolicUtils.to_poly!(poly_to_bs, bs_to_poly, resid / den)
                if newpoly isa SymbolicUtils.PolyVarT
                    newpoly = MP.polynomial(newpoly, SymbolicUtils.PolyCoeffT)
                end
                for (monomial, coeff) in result
                    if den_has_vars
                        monomial = monomial / den
                        @match monomial begin
                            BSImpl.Div(; num = n2, den = d2) => begin
                                den_has_vars2 = isequal(den, d2) || SymbolicUtils.query(in(vars), den)
                                if den_has_vars2
                                    coeff *= monomial
                                    monomial = 1
                                end
                            end
                            _ => nothing
                        end
                    else
                        coeff = coeff / den
                    end
                    mono = SymbolicUtils.to_poly!(poly_to_bs, bs_to_poly, monomial, false)
                    coeff = if SymbolicUtils.isconst(coeff)
                        unwrap_const(coeff)
                    else
                        SymbolicUtils.basicsymbolic_to_polyvar(bs_to_poly, coeff)
                    end
                    MA.operate!(+, newpoly, mono * coeff)
                end
                push!(subskeys, pvar)
                push!(subsvals, newpoly)
            end
            _ => nothing
        end
    end
    return MP.subs(poly, subskeys => subsvals)
end

function _semipoly_query_recurse(ex::SymbolicT)
    @match ex begin
        BSImpl.Term(; f, args) && if f === getindex && hasmetadata(args[1], VariableSource) end => false
        _ => true
    end
end

"""
$(TYPEDSIGNATURES)

Returns a tuple of two objects:

1. A dictionary of coefficients keyed by monomials in `vars` upto the given `degree`,
2. A residual expression which has all terms not represented as a product of monomial and a coefficient

`degree` should be a nonnegative number.

If  `consts` is set to `true`, then the returned dictionary will contain
a key `1` and the corresponding value will be the constant term. If `false`, the constant term will be part of the residual.
"""
function semipolynomial_form(expr, vars, degree::Real; consts = true)
    if degree < 0
        @warn "Degree for semi-polynomial form should be ≥ 0"
        return SemipolyDictT(), expr
    end
    vars = Set([unwrap(x) for x in vars])
    for v in vars
        v isa BasicSymbolic{VartypeT} || continue
        @match v begin
            BSImpl.Term(; f, args) && if f === getindex end => push!(vars, args[1])
            _ => nothing 
        end
    end
    expr = unwrap(expr)
    expr = expand(expr, false)
    poly_to_bs = Bijections.Bijection{SymbolicUtils.PolyVarT, BasicSymbolic{VartypeT}}()
    bs_to_poly = Bijections.active_inv(poly_to_bs)
    poly = SymbolicUtils.to_poly!(poly_to_bs, bs_to_poly, expr, false)
    poly = canonicalize_poly(poly_to_bs, bs_to_poly, poly, degree, vars)
    pvars = MP.variables(poly)
    nonpoly_mask = falses(length(pvars))
    in_vars_mask = falses(length(pvars))
    for (i, pvar) in enumerate(pvars)
        var = poly_to_bs[pvar]
        in_vars_mask[i] = var in vars
        in_vars_mask[i] && continue
        nonpoly_mask[i] = SymbolicUtils.query(in(vars), var; recurse = _semipoly_query_recurse)
    end
    result = SemipolyDictT()
    constant = SymbolicUtils.zeropoly()
    residual = SymbolicUtils.zeropoly()
    for t in MP.terms(poly)
        is_monomial_in_vars = true
        monomial_degree = 0
        true_monomial = SymbolicUtils.MonomialT()
        coeff_monomial = SymbolicUtils.MonomialT()
        for (i, exp) in enumerate(MP.exponents(t))
            iszero(exp) && continue
            monomial_degree += in_vars_mask[i] * exp
            MA.operate!(*, in_vars_mask[i] ? true_monomial : coeff_monomial, MP.variables(t)[i] ^ exp)
            is_monomial_in_vars &= !nonpoly_mask[i]
        end
        if is_monomial_in_vars && monomial_degree <= degree
            if monomial_degree == 0
                MA.operate!(+, constant, t)
            else
                mono = SymbolicUtils.from_poly(poly_to_bs, true_monomial)
                result[mono] = get(result, mono, 0) + MP.coefficient(t) * SymbolicUtils.from_poly(poly_to_bs, coeff_monomial)
            end
        else
            MA.operate!(+, residual, t)
        end
    end

    if consts && !iszero(constant)
        result[Const{VartypeT}(1)] = SymbolicUtils.from_poly(poly_to_bs, constant)
    else
        MA.operate!(+, residual, constant)
    end
    residual = SymbolicUtils.from_poly(poly_to_bs, residual)
    return result, residual
end

"""
$(TYPEDSIGNATURES)

For every expression in `exprs` computes the semi-polynomial form and
returns a tuple of two objects -- a vector of coefficient dictionaries,
and a vector of residual terms.

If  `consts` is set to `true`, then the returned dictionary will contain
a key `1` and the corresponding value will be the constant term. If `false`, the constant term will be part of the residual.
"""
function semipolynomial_form(exprs::AbstractArray, vars, degree::Real; consts = true)
    if degree < 0
        @warn "Degree for semi-polynomial form should be ≥ 0"
        return fill(SemipolyDictT(), size(exprs)), exprs
    end
    if any(iswrapped, vars)
        vars = map(unwrap, vars)
    end
    if !(vars isa AbstractSet)
        vars = Set(vars)
    end
    results = similar(exprs, SemipolyDictT)
    residuals = similar(exprs, BasicSymbolic{VartypeT})

    for (i, expr) in enumerate(exprs)
        results[i], residuals[i] = semipolynomial_form(expr, vars, degree; consts)
    end
    return results, residuals
end

"""
$(SIGNATURES)

Find coefficients of a polynomial in `vars`.

Returns a tuple of two elements:
1. A dictionary of coefficients keyed by monomials in `vars`
2. A residual expression which is the constant term

(Same as `semipolynomial_form(expr, vars, Inf)`)
"""
polynomial_coeffs(expr, vars) = semipolynomial_form(expr, vars, Inf)

"""
$(TYPEDSIGNATURES)

Returns a tuple of a sparse matrix `A`, and a residual vector `c` such that,
`A * vars + c` is the same as `exprs`.
"""
function semilinear_form(exprs::AbstractArray, vars)
    ds, nls = semipolynomial_form(exprs, vars, 1; consts = false)

    idxmap = Dict(v=>i for (i, v) in enumerate(vars))

    I = Int[]
    J = Int[]
    V = Num[]

    for (i, d) in enumerate(ds)
        for (k, v) in d
            push!(I, i)
            push!(J, idxmap[k])
            push!(V, v)
        end
    end

    sparse(I,J,V, length(exprs), length(vars)), wrap.(nls)
end

"""
$(TYPEDSIGNATURES)

Returns a tuple of 4 objects:

1. a matrix `A` of dimensions (m x n)
2. a matrix `B` of dimensions (m x (n+1)*n/2)
3. a vector `v2` of length (n+1)*n/2 containing monomials of `vars` upto degree 2 and zero where they are not required.
4. a residual vector `c` of length m.

where `n == length(exprs)` and `m == length(vars)`.


The result is arranged such that, `A * vars + B * v2 + c` is the same as `exprs`.
"""
function semiquadratic_form(exprs, vars)
    ds, nls = semipolynomial_form(exprs, vars, 2; consts = false)

    idxmap = Dict(v=>i for (i, v) in enumerate(vars))

    m, n = length(exprs), length(vars)
    I1 = Int[]
    J1 = Int[]
    V1 = Num[]

    I2 = Int[]
    J2 = Int[]
    V2 = Num[]

    v2_I = Int[]
    v2_V = Num[]

    for (i, d) in enumerate(ds)
        for (k, v) in d
            if pdegree(k) == 1
                push!(I1, i)
                push!(J1, idxmap[k])
                push!(V1, v)
            elseif pdegree(k) == 2
                push!(I2, i)
                if isop(k, ^)
                    b, e = arguments(k)
                    e = unwrap_const(e)
                    @assert e == 2
                    q = idxmap[b]
                    j = div(q*(q+1), 2)
                    push!(J2, j) # or div(q*(q-1), 2) + q
                    push!(V2, v)
                else
                    @assert isop(k, *)
                    a, b = arguments(k)
                    p, q = extrema((idxmap[a], idxmap[b]))
                    j = div(q*(q-1), 2) + p
                    push!(J2, j)
                    push!(V2, v)
                end
                push!(v2_I, j)
                push!(v2_V, k)
            else
                error("This should never happen")
            end
        end
    end


    #v2 = SparseVector(div(n * (n + 1), 2), v2_I, v2_V) # When it works in the future
    # until then
    v2 = zeros(Num, div(n * (n + 1), 2))
    v2[v2_I] .= v2_V

    tuple(sparse(I1,J1,V1, m, n),
          sparse(I2,J2,V2, m, div(n * (n + 1), 2)),
          v2,
          wrap.(nls))
end

isop(x, op) = iscall(x) && operation(x) === op
isop(op) = Base.Fix2(isop, op)

function pdegrees(x)
    x = unwrap(x)
    if ismul(x)
        return x.dict
    elseif isdiv(x)
        num_dict = pdegrees(x.num)
        den_dict = pdegrees(x.den)
        inv_den_dict = Dict(keys(den_dict) .=> map(-, values(den_dict)))
        mergewith(+, num_dict, inv_den_dict)
    elseif ispow(x)
        base, exp = arguments(x)
        dict = pdegrees(base)
        degrees = map(degree -> degree * unwrap_const(exp), values(dict))
        Dict(keys(dict) .=> degrees)
    elseif issym(x) || iscall(x)
        return Dict(x=>1)
    elseif SymbolicUtils.isconst(x) || x isa Number
        return Dict()
    else
        error("pdegrees for $x unknown")
    end
end

pdegree(x::Number) = 0
function pdegree(x)
    degree_dict = pdegrees(x)
    if isempty(degree_dict)
        return 0
    end
    sum(values(degree_dict))
end
