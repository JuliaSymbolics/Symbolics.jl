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

This notation transforms multiplication into addition of exponent vertors, division into
subtraction, exponentiation into addition.

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
function tomonomial(m::SemiMonomial{T}, vars)::Union{Real, Symbolic} where {T}
    indices = findall(x -> x > 0, m.degrees)
    dict = Dict(vars[i] => Int(m.degrees[i]) for i in indices)
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

function init_semipoly_vars(vars)
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
    semilinear_form(exprs::AbstractVector, vars::AbstractVector, [consts=false])

Returns a tuple of a sparse matrix `A`, and a residual vector `c` such that,

`A * vars + c` is the same as `exprs`.
If `consts` = true then, `A * [1, vars] + c == exprs`.
Here `c` will only contain nonlinear terms and not constant terms.
"""
function semilinear_form(exprs::AbstractArray, vars, consts=false)
    exprs = unwrap.(exprs)
    vars = init_semipoly_vars(vars)
    ds, nls = semipolynomial_form(exprs, vars, 1, consts)

    idxmap = Dict(v=>i for (i, v) in enumerate(consts ? [1, vars...] : vars))

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

    sparse(I,J,V, length(exprs), length(vars) + consts), wrap.(nls)
end

"""
    semiquadratic_form(exprs::AbstractVector, vars::AbstractVector)

Returns a tuple of 4 objects:

1. a matrix `A` of dimensions (m x n)
2. a matrix `B` of dimensions (m x (n+1)*n/2)
3. a vector `v2` of length (n+1)*n/2 containing monomials of `vars` upto degree 2 and zero where they are not required.
4. a residual vector `c` of length m.

where `n == length(exprs)` and `m == length(vars)`.


The result is arranged such that, `A * vars + B * v2 + c` is the same as `exprs`.
"""
function semiquadratic_form(exprs, vars)
    exprs = unwrap.(exprs)
    vars = init_semipoly_vars(vars)
    ds, nls = semipolynomial_form(exprs, vars, 2)

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
                    @assert e == 2
                    q = idxmap[b]
                    j = div(q*(q+1), 2)
                    push!(J2, j) # or div(q*(q-1), 2) + q
                    push!(V2, v)
                else
                    @assert isop(k, *)
                    a, b = unsorted_arguments(k)
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


## Utilities

function partition(n, parts)
    parts == 1 && return n
    [[[i, p...] for p in partition(n-i, parts-1)]
     for i=0:n] |> Iterators.flatten |> collect
end

unwrap_bp(x::BoundedDegreeMonomial) = x.p * x.coeff
function unwrap_bp(x)
    x = unwrap(x)
    istree(x) ? similarterm(x, operation(x), map(unwrap_bp, unsorted_arguments(x))) : x
end

cautious_sum(nls) = isempty(nls) ? 0 : isone(length(nls)) ? unwrap_bp(first(nls)) : sum(unwrap_bp, nls)

_mul(x) = isempty(x) ? 1 : isone(length(x)) ? first(x) : prod(x)

function mul(a, b, deg)
    if isop(a, +)
        return Term{symtype(a)}(+, map(x->mul(x, b, deg), unsorted_arguments(a)))
    elseif isop(b, +)
        return Term{symtype(a)}(+, mul.((a,), unsorted_arguments(b), deg))
    elseif a isa BoundedDegreeMonomial
        return BoundedDegreeMonomial(a.p, a.coeff * b, a.overdegree)
    elseif b isa BoundedDegreeMonomial
        return BoundedDegreeMonomial(b.p, a * b.coeff, b.overdegree)
    else
        return a * b
    end
end

function mul(a::BoundedDegreeMonomial, b::BoundedDegreeMonomial, deg)
    if a.overdegree || b.overdegree || pdegree(a.p) + pdegree(b.p) > deg
        highdegree(a.p * b.p * a.coeff * b.coeff)
    else
        BoundedDegreeMonomial(a.p * b.p, a.coeff * b.coeff, false)
    end
end

function mul_bounded(xs, deg)
    length(xs) == 1 ?
        first(xs) :
        reduce((x,y) -> mul(x,y, deg), xs)
end

function pow(a::BoundedDegreeMonomial, b, deg)
    return BoundedDegreeMonomial((a.p)^b, a.coeff ^ b, !_isone(a.p) && b > deg)
end

pow(a, b, deg) = a^b

function pow_of_add(a, b, deg, vars, consts)
    b == 0 && return 1
    b == 1 && return a

    @assert isop(a, +)
    within_deg(x) = (d=_degree(a);d <= deg)

    a = mark_and_exponentiate(a, vars, deg, consts)
    args = all_terms(a)
    mindeg = minimum(_degree, args)
    maxdeg = maximum(_degree, args)

   #if maxdeg * b <= deg
   #    return mark_and_exponentiate(expand(unwrap_bp(a^b)), vars, deg)
   #end

    if mindeg * b > deg
        return highdegree(unwrap_bp(a^b)) # That's pretty high degree
    end

    interesting = filter(within_deg, args)
    nls = filter(!within_deg, args)

    if isempty(interesting)
        a^b # Don't need to do anything!
    else
        q = partial_multinomial_expansion(interesting, b, deg, consts)
        if maxdeg * b <= deg
            q # q is the whole enchilada
        else
            Term{Real}(+, [all_terms(q)...,
                           map(highdegree, all_terms(unwrap_bp(a^b - q)))...])
        end
    end
end

function partial_multinomial_expansion(xs, exp, deg, consts)
    zs = filter(iszero∘_degree, xs)
    nzs = filter((!iszero)∘_degree, xs)
    function add_consts(xs)
        !consts && return +(xs...)
        s = 0
        for x in xs
            if x isa BoundedDegreeMonomial
                @assert x.p == 1
                s += x.coeff
            else
                s += x
            end
        end
        return BoundedDegreeMonomial(1, s, false)
    end
    if isempty(zs)
        terms = nzs
    else
        terms = [add_consts(zs), nzs...]
    end
    degs = map(_degree, terms)


    q = []
    nfact = factorial(exp)
    for ks in partition(exp, length(terms))
        td = sum(degs .* ks)

        coeff = div(nfact, prod(factorial, ks))
        function monomial(coeff, terms, ks)
            mul_bounded([coeff, (pow(x,y, deg) for (x, y)
                                 in zip(terms, ks)
                                 if !iszero(y))...], Inf)
        end
        if td == 0
            m = monomial(coeff, terms, ks)
            push!(q, monomial(coeff, terms, ks))
        elseif td <= deg
            push!(q,
                  mul_bounded(vcat(coeff,
                                   [pow(x,y, deg)
                                    for (x, y) in zip(terms, ks) if !iszero(y)]),
                              deg))
        end
    end
    return Term{Real}(+, q)
end

# used to get all arguments of a possibly nested `Term` with + operation or `Add`.
function all_terms(x)
    if istree(x) && operation(x) == (+)
        collect(Iterators.flatten(map(all_terms, unsorted_arguments(x))))
    else
        (x,)
    end
end
