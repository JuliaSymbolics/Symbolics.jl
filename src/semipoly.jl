import SymbolicUtils.Rewriters: RestartedChain
using DataStructures

export semipolynomial_form, semilinear_form, semiquadratic_form, polynomial_coeffs

"""
$(TYPEDEF)

A compact notation for monomials and also non-monomial terms.

# Attrtibutes
$(TYPEDFIELDS)

For example, the monomial ``2 x^3 y^5 z^7`` about the variables ``(x, y, z)`` is simply
expressed as `SemiMonomial(coeff = 2, degrees = [3, 5, 7])`.

This struct is called *semi* because it can also represent non-monomial terms.
For example, ``5 b^{2.5} \\tan(c) / a^{\\frac12}`` about ``(a, b)`` is
`SemiMonomial(coeff = 5tan(c), degrees = [-1//2, 2.5])`.
Note that here ``c`` is treated as a constant.

This notation transforms multiplication of monomials into the addition of exponent vertors,
and division into subtraction. The exponentiation with a monomial as the base and a real
number as the exponent is tranformed to the multiplication of the exponent vector and a
scalar.

The parametric type `T` depends on the types of the associated variables. For example,
when the variables are declared using `@variables x::Int32 y::Int64 z::Rational{Int32}`,
`T` should be `Rational{Int64}` derived by `promote_type(Symbolics.symtype.([x, y, z])...)`.

See also
[Wikipedia: Monomial - Multi-index notation](https://en.wikipedia.org/wiki/Monomial#Multi-index_notation).
"""
struct SemiMonomial{T}
    "coefficient"
    coeff::Any
    "exponent vector"
    degrees::Vector{N} where {N <: Real}
end

Base.:+(a::SemiMonomial) = a
function Base.:+(a::SemiMonomial{S}, b::SemiMonomial{T}) where {S, T}
    SymbolicUtils.Term{promote_symtype(+, S, T)}(+, [a, b])
end
function Base.:+(a::SymbolicUtils.Term{S, M}, b::SemiMonomial{T}) where {S, T, M}
    SymbolicUtils.Term{promote_symtype(+, S, T)}(+, [a.arguments; b])
end
Base.:+(a::SemiMonomial{T}, b::SymbolicUtils.Term{S, M}) where {S, T, M} = b + a

Base.:*(m::SemiMonomial) = m
function Base.:*(a::SemiMonomial{S}, b::SemiMonomial{T}) where {S, T}
    SemiMonomial{promote_symtype(*, S, T)}(a.coeff * b.coeff, a.degrees + b.degrees)
end
function Base.:*(m::SemiMonomial{T}, t) where {T}
    if istree(t) && operation(t) == (+)
        args = collect(all_terms(t))
        return SymbolicUtils.Term(+, (m,) .* args)
    end
    SemiMonomial{promote_symtype(*, T, symtype(t))}(m.coeff * t, m.degrees)
end
Base.:*(t, m::SemiMonomial) = m * t

function Base.:/(a::SemiMonomial{S}, b::SemiMonomial{T}) where {S, T}
    SemiMonomial{promote_symtype(/, S, T)}(a.coeff / b.coeff, a.degrees - b.degrees)
end

function Base.:^(base::SemiMonomial{T}, exp::Real) where {T}
    SemiMonomial{promote_symtype(^, T, typeof(exp))}(base.coeff^exp, base.degrees * exp)
end

"""
$(SIGNATURES)

Check if `x` is of type [`SemiMonomial`](@ref).
"""
issemimonomial(x) = x isa SemiMonomial

"""
$(SIGNATURES)

Return true if `m` is a [`SemiMonomial`](@ref) and satisfies the definition of a monomial.

A monomial, also called power product, is a product of powers of variables with nonnegative
integer exponents.

See also [Wikipedia: Monomial](https://en.wikipedia.org/wiki/Monomial).
"""
function ismonomial(m, vars)::Bool
    if !(m isa SemiMonomial)
        return false
    end
    for degree in m.degrees
        if !isinteger(degree) || degree < 0
            return false
        end
    end
    !has_vars(m.coeff, vars)
end

"""
$(TYPEDSIGNATURES)

Return true is `m` is a [`SemiMonomial`](@ref), satisfies the definition of a monomial and
its degree is less than or equal to `degree_bound`. See also [`ismonomial`](@ref).
"""
function isboundedmonomial(m, vars, degree_bound::Real)::Bool
    ismonomial(m, vars) && _degree(m) <= degree_bound
end

"""
$(SIGNATURES)

Construct a [`SemiMonomial`](@ref) object with `expr` as its coefficient and 0 degrees.
"""
function non_monomial(expr, vars)::SemiMonomial
    degrees = zeros(Int, length(vars))
    SemiMonomial{symtype(expr)}(expr, degrees)
end

"""
$(TYPEDSIGNATURES)

Return true if the degrees of `m` are all 0s and its coefficient is a `Real`.
"""
function Base.:isreal(m::SemiMonomial)::Bool
    _degree(m) == 0 && unwrap(m.coeff) isa Real
end
Base.:isreal(::Symbolic) = false

"""
$(TYPEDSIGNATURES)

Transform `m` to a `Real`.

Assume `isreal(m) == true`, otherwise calling this function does not make sense.
"""
function Base.:real(m::SemiMonomial)::Real
    if isinteger(m.coeff)
        return Int(m.coeff)
    end
    return m.coeff
end

# needed for `SymbolicUtils.expand`
symtype(::SemiMonomial{T}) where {T} = T

TermInterface.issym(::SemiMonomial) = true

Base.:nameof(m::SemiMonomial) = Symbol(:SemiMonomial, m.coeff, m.degrees)

isop(x, op) = istree(x) && operation(x) === op
isop(op) = Base.Fix2(isop, op)

pdegree(x::Mul) = sum(values(x.dict))
pdegree(x::Union{Sym, Term}) = 1
pdegree(x::Pow) = pdegree(x.base) * x.exp
pdegree(x::Number) = 0

_degree(x::SemiMonomial) = sum(x.degrees)
_degree(x::Symbolic) = isop(x, *) ? sum(_degree, unsorted_arguments(x)) : 0

bareterm(x, f, args; kw...) = Term{symtype(x)}(f, args)

function mark_and_exponentiate(expr, vars)
    # Step 1
    # Mark all the interesting variables -- substitute without recursing into nl forms
    expr′ = mark_vars(expr, vars)

    # Step 2
    # Construct and propagate BoundedDegreeMonomial for ^ and * and /

    # does not do fraction simplification
    rules = [@rule (~a::issemimonomial)^(~b::isreal) => (~a)^real(~b)
             @rule (~a::isop(+))^(~b::isreal) => expand((~a)^real(~b))
             @rule *(~~xs) => expand(*(~~xs...))
             @rule (~a::issemimonomial) / (~b::issemimonomial) => (~a) / (~b)]
    expr′ = Postwalk(RestartedChain(rules), similarterm = bareterm)(expr′)
end

function semipolyform_terms(expr, vars)
    expr = mark_and_exponentiate(expr, vars)
    if istree(expr) && operation(expr) == (+)
        args = collect(all_terms(expr))
        return args
    elseif isreal(expr) && iszero(real(expr)) # when `expr` is just a 0
        return []
    else
        return [expr]
    end
end
semipolyform_terms(vars) = Base.Fix2(semipolyform_terms, vars)

"""
$(TYPEDSIGNATURES)

Return true if `expr` contains any variables in `vars`.
"""
function has_vars(expr, vars)::Bool
    if expr in vars
        return true
    elseif istree(expr)
        for arg in unsorted_arguments(expr)
            if has_vars(arg, vars)
                return true
            end
        end
    end
    return false
end

function mark_vars(expr, vars)
    index = findfirst(isequal(expr), vars)
    if !isnothing(index)
        degrees = zeros(Int, length(vars))
        degrees[index] = 1
        return SemiMonomial{symtype(expr)}(1, degrees)
    elseif !istree(expr)
        return non_monomial(expr, vars)
    end
    op = operation(expr)
    args = arguments(expr)
    if op === (+) || op === (*)
        return Term{symtype(expr)}(op, map(mark_vars(vars), args))
    elseif op === (^) || op == (/)
        @assert length(args) == 2
        return Term{symtype(expr)}(op, map(mark_vars(vars), args))
    elseif length(args) == 1
        if op == sqrt
            base = mark_vars(args[1], vars)
            degrees = zeros(Int, length(vars))
            exp = SemiMonomial{Rational{Int}}(1 // 2, degrees)
            return Pow(base, exp)
        elseif linearity_1(op)
            return Term{symtype(expr)}(op, mark_vars(args[1], vars))
        end
    end
    return non_monomial(expr, vars)
end
mark_vars(vars) = Base.Fix2(mark_vars, vars)

"""
$(TYPEDSIGNATURES)

Transform `m` into its corresponding `Real` number or `SymbolicUtils.Symbolic` form.

Assume `m` satisfies [`ismonomial`](@ref).
"""
function tomonomial(m::SemiMonomial{T}, vars::OrderedSet)::Union{Real, Symbolic} where {T}
    indices = findall(x -> x > 0, m.degrees)
    dict = Dict(vars.dict.keys[i] => Int(m.degrees[i]) for i in indices)
    Mul(T, 1, dict)
end

# Transform `SemiMonomial` and `SymbolicUtils.Symbolic` to their corresponding appropriate
# `SymbolicUtils.Symbolic` subtypes.
function unwrap_sm(m::SemiMonomial{T}, vars) where {T}
    dict_positive = Dict()
    sizehint!(dict_positive, length(vars))
    # deal with negative powers separately to avoid making negative degrees
    # in `Mul` or `Pow`
    dict_negative = Dict()
    sizehint!(dict_positive, length(vars))
    for (var, degree) in zip(vars, m.degrees)
        if isinteger(degree)
            degree = Int(degree)
        end
        if degree > 0
            dict_positive[var] = degree
        elseif degree < 0
            dict_negative[var] = -degree
        end
    end
    positive = if m.coeff isa Number
        Mul(T, m.coeff, dict_positive)
    else
        m.coeff * Mul(T, 1, dict_positive)
    end
    negative = Mul(T, 1, dict_negative)
    positive / negative
end
function unwrap_sm(x, vars)
    x = unwrap(x)
    if istree(x)
        similarterm(x, operation(x), map(unwrap_sm(vars), arguments(x)))
    else
        x
    end
end
unwrap_sm(vars) = Base.Fix2(unwrap_sm, vars)

function bifurcate_terms(terms, vars, degree_bound::Real)
    # Step 4: Bifurcate polynomial and nonlinear parts:

    monomials = filter(arg -> isboundedmonomial(arg, vars, degree_bound), terms)

    polys_dict = Dict()
    sizehint!(polys_dict, length(monomials))
    for m in monomials
        monomial = tomonomial(m, vars)
        if haskey(polys_dict, monomial)
            polys_dict[monomial] += m.coeff
        else
            polys_dict[monomial] = m.coeff
        end
    end
    if length(monomials) == length(terms)
        return polys_dict, 0
    end
    nl_terms = setdiff(terms, monomials)
    nl = unwrap_sm(sum(unwrap_sm(vars), nl_terms), vars)
    return polys_dict, nl
end

function init_semipoly_vars(vars)::OrderedSet
    set = OrderedSet(unwrap.(vars))
    @assert length(set) == length(vars) # vars passed to semi-polynomial form must be unique
    set
end

"""
$(TYPEDSIGNATURES)

Returns a tuple of two objects:

1. A dictionary of coefficients keyed by monomials in `vars` upto the given `degree`,
2. A residual expression which has all terms not represented as a product of monomial and a coefficient

`degree` should be a nonnegative number.

See also
[Wikipedia: Polynomial](https://en.wikipedia.org/wiki/Polynomial).
"""
function semipolynomial_form(expr, vars, degree::Real)
    if degree < 0
        @warn "Degree for semi-polynomial form should be ≥ 0"
        return Dict(), expr
    end
    vars = init_semipoly_vars(vars)
    expr = unwrap(expr)
    terms = semipolyform_terms(expr, vars)
    bifurcate_terms(terms, vars, degree)
end

"""
$(TYPEDSIGNATURES)

For every expression in `exprs` computes the semi-polynomial form and
returns a tuple of two objects -- a vector of coefficient dictionaries,
and a vector of residual terms.
"""
function semipolynomial_form(exprs::AbstractArray, vars, degree::Real)
    if degree < 0
        @warn "Degree for semi-polynomial form should be ≥ 0"
        return fill(Dict(), length), exprs
    end
    vars = init_semipoly_vars(vars)
    exprs = unwrap.(exprs)
    matches = map(semipolyform_terms(vars), exprs)
    tmp = map(match -> bifurcate_terms(match, vars, degree), matches)
    map(first, tmp), map(last, tmp)
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
function semilinear_form(exprs::AbstractArray, vars::AbstractVector)
    vars = init_semipoly_vars(vars)
    exprs = unwrap.(exprs)

    matches = map(semipolyform_terms(vars), exprs)

    I = Int[]
    J = Int[]
    V = Num[]

    nls = Vector{Union{Real, SymbolicUtils.Symbolic}}(undef, length(exprs))
    nls .= 0
    for (row_index, terms) in enumerate(matches)
        constant_linear_terms = filter(t -> isboundedmonomial(t, vars, 1), terms)
        for term in constant_linear_terms
            if _degree(term) == 1 # linear term
                col_index = findfirst(Bool.(term.degrees))
                push!(I, row_index)
                push!(J, col_index)
                push!(V, term.coeff)
            else # constant term
                nls[row_index] += unwrap_sm(term, vars)
            end
        end
        nonlinear_terms = setdiff(terms, constant_linear_terms)
        if !isempty(nonlinear_terms)
            nls[row_index] += sum(unwrap_sm(vars), nonlinear_terms)
        end
    end

    sparse(I, J, V, length(exprs), length(vars)), wrap.(nls)
end

"""
$(TYPEDSIGNATURES)

Returns a tuple of 4 objects:

1. a matrix `A` of dimensions (m x n)
2. a matrix `B` of dimensions (m x (n+1)*n/2)
3. a vector `v2` of length (n+1)*n/2 containing monomials of `vars` upto degree 2 and zero where they are not required.
4. a residual vector `c` of length m.

where `m == length(exprs)` and `n == length(vars)`.

The result is arranged such that, `A * vars + B * v2 + c` is the same as `exprs`.
"""
function semiquadratic_form(exprs::AbstractVector, vars::AbstractVector)
    vars = init_semipoly_vars(vars)
    exprs = unwrap.(exprs)

    matches = map(semipolyform_terms(vars), exprs)

    I1 = Int[]
    J1 = Int[]
    V1 = Num[]

    I2 = Int[]
    J2 = Int[]
    V2 = Num[]

    v2_I = Int[]
    v2_V = Num[]

    nls = Vector{Union{Real, SymbolicUtils.Symbolic}}(undef, length(exprs))
    nls .= 0
    for (row_index, terms) in enumerate(matches)
        const_linear_quadratic_terms = filter(t -> isboundedmonomial(t, vars, 2), terms)
        for term in const_linear_quadratic_terms
            degree = _degree(term)
            if degree == 1 # linear term
                col_index = findfirst(Bool.(term.degrees))
                push!(I1, row_index)
                push!(J1, col_index)
                push!(V1, term.coeff)
            elseif degree == 2 # quadratic term
                push!(I2, row_index)
                push!(V2, term.coeff)
                var_indices = findall(Bool.(sign.(term.degrees)))
                if length(var_indices) == 1 # of the form x²
                    var_index = var_indices[1]
                    col_index = var_index * (var_index + 1) ÷ 2
                    push!(v2_V, vars.dict.keys[var_index]^2)
                else # of the form x * y
                    index₁, index₂ = extrema(var_indices)
                    col_index = index₂ * (index₂ - 1) ÷ 2 + index₁
                    push!(v2_V, vars.dict.keys[index₁] * vars.dict.keys[index₂])
                end
                push!(J2, col_index)
                push!(v2_I, col_index)
            else # constant term
                nls[row_index] += unwrap_sm(term, vars)
            end
        end
        # nonquadratic_terms = filter(t -> !isboundedmonomial(t, vars, 2), terms)
        nonquadratic_terms = setdiff(terms, const_linear_quadratic_terms)
        if !isempty(nonquadratic_terms)
            nls[row_index] += sum(unwrap_sm(vars), nonquadratic_terms)
        end
    end
    m = length(exprs)
    n = length(vars)

    #v2 = SparseVector(div(n * (n + 1), 2), v2_I, v2_V) # When it works in the future
    # until then
    v2 = zeros(Num, div(n * (n + 1), 2))
    v2[v2_I] .= v2_V

    tuple(sparse(I1, J1, V1, m, n),
          sparse(I2, J2, V2, m, div(n * (n + 1), 2)),
          v2,
          wrap.(nls))
end

# used to get all arguments of a possibly nested `Term` with + operation or `Add`.
function all_terms(x)
    if istree(x) && operation(x) == (+)
        collect(Iterators.flatten(map(all_terms, unsorted_arguments(x))))
    else
        (x,)
    end
end
