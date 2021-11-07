import SymbolicUtils: @ordered_acrule, unpolyize
using DataStructures

struct BoundedDegreeMonomial
    p::Union{Mul, Pow, Int, Sym}
    coeff::Any
    overdegree::Bool
end
SymbolicUtils.symtype(b::BoundedDegreeMonomial) = symtype(b.p)

# required by unsorted_arguments etc.
SymbolicUtils.unstable_pow(x::BoundedDegreeMonomial, i::Integer) = (@assert(i==1); x)


function flatten_upto_deg(expr, deg)
    isbp(x) = x isa BoundedDegreeMonomial && !x.overdegree
    st(x, f, args;kw...) = Term{symtype(x)}(f, args)
    rules = [@rule (~a::isbp) ^ (~b::(x-> x isa Integer)) =>
             BoundedDegreeMonomial(((~a).p)^(~b), (~a).coeff ^ (~b), ~b > deg)
             @rule *(~~x) => mul_bounded(~~x, deg)]
    Postwalk(RestartedChain(rules), similarterm=st)(expr)
end

function semipolyform_terms(expr, vars::OrderedSet, deg)
    # Step 1
    # Mark all the interesting variables -- substitute without recursing into nl forms

    expr′ = mark_vars(expr, vars)


    # Step 2
    # Construct and propagate BoundedDegreeMonomial for ^ and *

    isbp(x) = x isa BoundedDegreeMonomial && !x.overdegree

    st(x, f, args;kw...) = Term{symtype(x)}(f, args)


    expr′ = flatten_upto_deg(expr′, deg)

    # Step 3: every term now has at most one valid monomial -- find the coefficient of it.
    expr′ = Postwalk(Chain([@rule *(~~x, ~a::isbp, ~~y) =>
                                BoundedDegreeMonomial((~a).p,
                                                      (~a).coeff * _mul(~~x) * _mul(~~y),
                                                      (~a).overdegree)]),
                     similarterm=st)(expr′)
end

all_terms(x) = istree(x) && operation(x) == (+) ? Iterators.flatten(map(all_terms, unsorted_arguments(x))) : (x,)

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

unwrap_bp(x::BoundedDegreeMonomial) = x.p * x.coeff
unwrap_bp(x) = x

cautious_sum(nls) = isempty(nls) ? 0 : isone(length(nls)) ? unwrap_bp(first(nls)) : sum(unwrap_bp, nls)

_mul(x) = isempty(x) ? 1 : isone(length(x)) ? first(x) : prod(x)

function mul_bounded(xs, deg)
    length(xs) == 1 ?
        first(xs) :
        reduce((x,y) -> mul(x,y, deg), xs)
end

pdegree(x::Mul) = sum(values(x.dict))
pdegree(x::Sym) = 1
pdegree(x::Pow) = pdegree(x.base) * x.exp
pdegree(x::Number) = 0

_deg(x) = pdegree(x)
_deg(x::Add) = maximum(_deg, unsorted_arguments(x))

function mul(a, b, deg)
    isop(x, op) = istree(x) && operation(x) === op
    if isop(a, +)
        return Term{symtype(a)}(+, mul.(unsorted_arguments(a), (b,), deg))
    elseif isop(b, +)
        return Term{symtype(a)}(+, mul.((a,), unsorted_arguments(b), deg))
    elseif a isa BoundedDegreeMonomial
        return BoundedDegreeMonomial(a.p, a.coeff * b, false)
    elseif b isa BoundedDegreeMonomial
        return BoundedDegreeMonomial(b.p, a * b.coeff, false)
    else
        return a * b
    end
end

function mul(a::BoundedDegreeMonomial, b::BoundedDegreeMonomial, deg)
    if a.overdegree || b.overdegree || pdegree(a.p) + pdegree(b.p) > deg
        BoundedDegreeMonomial(1, a.p * b.p * a.coeff * b.coeff, true)
    else
        BoundedDegreeMonomial(a.p * b.p, a.coeff * b.coeff, false)
    end
end


function mark_vars(expr, vars)
    if !istree(expr)
        if expr in vars
            return BoundedDegreeMonomial(expr, 1, false)
        else
            return expr
        end
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
