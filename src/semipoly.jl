import SymbolicUtils.Rewriters: RestartedChain
using DataStructures

struct BoundedDegreeMonomial
    p::Union{Mul, Pow, Int, Sym, Term}
    coeff::Any
    overdegree::Bool
end

highdegree(x) = BoundedDegreeMonomial(1, x, true)
highdegree(x::BoundedDegreeMonomial) = (@assert(x.overdegree); x)
SymbolicUtils.symtype(b::BoundedDegreeMonomial) = symtype(b.p)

isop(x, op) = istree(x) && operation(x) === op
isop(op) = x -> isop(x, op)
# required by unsorted_arguments etc.
SymbolicUtils.unstable_pow(x::BoundedDegreeMonomial, i::Integer) = (@assert(i==1); x)

function semipolyform_terms(expr, vars::OrderedSet, deg)
    # Step 1
    # Mark all the interesting variables -- substitute without recursing into nl forms

    expr′ = mark_vars(expr, vars)


    # Step 2
    # Construct and propagate BoundedDegreeMonomial for ^ and *

    isbp(x) = x isa BoundedDegreeMonomial && !x.overdegree
    st(x, f, args;kw...) = Term{symtype(x)}(f, args)

    rules = [@rule (~a::isbp) ^ (~b::(x-> x isa Integer)) =>
             BoundedDegreeMonomial(((~a).p)^(~b), (~a).coeff ^ (~b), ~b > deg)
             @rule (~a::isop(+)) ^ (~b::(x -> x isa Integer)) => sum_pow(~a, ~b, deg)
             @rule *(~~x) => mul_bounded(~~x, deg)]

    expr′ = Postwalk(RestartedChain(rules), similarterm=st)(expr′)

    # Step 3: every term now has at most one valid monomial -- find the coefficient of it.
    expr′ = Postwalk(Chain([@rule *(~~x, ~a::isbp, ~~y) =>
                                BoundedDegreeMonomial((~a).p,
                                                      (~a).coeff * _mul(~~x) * _mul(~~y),
                                                      (~a).overdegree)]),
                     similarterm=st)(expr′)
end

_degree(x::BoundedDegreeMonomial) = x.overdegree ? Inf : pdegree(x.p)
_degree(x) = 0

function sum_pow(a, b, deg)
    b == 0 && return 1
    b == 1 && return a

    @assert isop(a, +)
    within_deg(x) = (d=_degree(a);d <= deg)

    args = all_terms(a)
    mindeg = minimum(_degree, args)
    maxdeg = maximum(_degree, args)

    if maxdeg * b <= deg
        return expand(a^b) # & TODO: this has BoundedDegreeMonomials so won't work, make it work nicely
    end

    if mindeg * b > deg
        return highdegree(a^b) # That's pretty high degree
    end

    interesting = filter(within_deg, args)
    nls = filter(!within_deg, args)

    if isempty(interesting)
        return a^b # Don't need to do anything!
    else
        q = partial_multinomial_expansion(interesting, b, deg)
        return Term{Real}(+, [q, a^b - unwrap_bp(q)])
    end
end

function partial_multinomial_expansion(xs, exp, deg)
    zs = filter(iszero∘_degree, xs)
    nzs = filter((!iszero)∘_degree, xs)
    if isempty(zs)
        terms = nzs
    else
        terms = [+(zs...), nzs...]
    end
    degs = map(_degree, terms)

    q = []
    for ks in partition(exp, length(terms))
        td = sum(degs .* ks)
        td == 0 && continue # We are not concerned about constants
        if td <= deg
            push!(q, mul_bounded([pow(x,y, deg) for (x, y) in zip(terms, ks) if !iszero(y)], deg))
        end
    end
    return Term{Real}(+, q)
end

function partition(n, parts)
    parts == 1 && return n
    [[[i, p...] for p in partition(n-i, parts-1)]
     for i=0:n] |> Iterators.flatten |> collect
end

all_terms(x) = istree(x) && operation(x) == (+) ? collect(Iterators.flatten(map(all_terms, unsorted_arguments(x)))) : (x,)

function bifurcate_terms(expr)
    # Step 4: Bifurcate polynomial and nonlinear parts:

    isbp(x) = x isa BoundedDegreeMonomial && !x.overdegree
    if isbp(expr)
        return Dict(expr.p => expr.coeff), 0
    elseif istree(expr) && operation(expr) == (+)
        args = collect(all_terms(expr))
        polys = filter(isbp, args)

        poly = Dict()
        for p in polys
            if haskey(poly, p.p)
                poly[p.p] += p.coeff
            else
                poly[p.p] = p.coeff
            end
        end

        nls = filter(!isbp, args)
        nl = cautious_sum(nls)

        return poly, nl
    else
        return Dict(), expr
    end
end

function init_semipoly_vars(vars)
    set = OrderedSet(unwrap.(vars))
    @assert length(set) == length(vars) "vars passed to semi-polynomial form must be unique"
    set
end

function semipolynomial_form(exprs::AbstractArray, vars, degree)
    exprs = unwrap.(exprs)
    vars = init_semipoly_vars(vars)

    matches = map(x -> semipolyform_terms(x, vars, degree), exprs)
    tmp = map(bifurcate_terms, matches)
    dicts, nls = map(first, tmp), map(last, tmp)
end

function semipolynomial_form(expr, vars, degree)
    expr = unwrap(expr)
    vars = init_semipoly_vars(vars)

    bifurcate_terms(semipolyform_terms(x, vars, degree))
end

function semilinear_form(exprs::AbstractArray, vars)
    exprs = unwrap.(exprs)
    vars = init_semipoly_vars(vars)
    ds, nls = semipolynomial_form(exprs, vars, 1)

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

unwrap_bp(x::BoundedDegreeMonomial) = x.p * x.coeff
unwrap_bp(x) = istree(x) ? similarterm(x,
                                       operation(x),
                                       map(unwrap_bp, unsorted_arguments(x))) : x

cautious_sum(nls) = isempty(nls) ? 0 : isone(length(nls)) ? unwrap_bp(first(nls)) : sum(unwrap_bp, nls)

_mul(x) = isempty(x) ? 1 : isone(length(x)) ? first(x) : prod(x)

function mul_bounded(xs, deg)
    length(xs) == 1 ?
        first(xs) :
        reduce((x,y) -> mul(x,y, deg), xs)
end

pdegree(x::Mul) = sum(values(x.dict))
pdegree(x::Union{Sym, Term}) = 1
pdegree(x::Pow) = pdegree(x.base) * x.exp
pdegree(x::Number) = 0

function mul(a, b, deg)
    if isop(a, +)
        return Term{symtype(a)}(+, mul.(unsorted_arguments(a), (b,), deg))
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

function pow(a::BoundedDegreeMonomial, b, deg)
    return BoundedDegreeMonomial((a.p)^b, a.coeff ^ b, b > deg)
end
pow(a, b, deg) = a^b


function mark_vars(expr, vars)

    if expr in vars
        return BoundedDegreeMonomial(expr, 1, false)
    end

    if !istree(expr)
        return expr
    else
        op = operation(expr)
        args = unsorted_arguments(expr)

        if op === (+) || op === (*)
            return Term{symtype(expr)}(op, map(x -> mark_vars(x, vars), args))
        elseif op === (^)
            base, exp = arguments(expr)
            Term{symtype(expr)}(^, [mark_vars(base, vars), exp])
        elseif length(args) == 1
            if linearity_1(op)
                return Term{symtype(expr)}(op, map(x -> mark_vars(x, vars), args))
            else
                return expr
            end
        else
            return expr
        end
    end
end
