import SymbolicUtils: @ordered_acrule, unpolyize
import DynamicPolynomials

struct BoundedDegreeMonomial
    p::Union{Mul, Pow, Int, Sym}
    coeff::Any
    overdegree::Bool
end

SymbolicUtils.unstable_pow(x::BoundedDegreeMonomial, i::Integer) = (@assert(i==1); x)

function Base.isequal(a::BoundedDegreeMonomial, b::BoundedDegreeMonomial)
    a === b && return true
    a.p === b.p && a.coeff === b.coeff && return true
    isequal(a.p, b.p) && isequal(a.coeff, b.coeff)
end

function semipolyform(expr, vars, deg)
    # Step 1
    # Mark all the interesting variables -- substitute without recursing into nl forms
    expr′ = mark_vars(expr, vars)


    # Step 2
    # Construct and propagate BoundedDegreeMonomial for ^ and *

    isbp(x) = x isa BoundedDegreeMonomial && !x.overdegree

    st(x, f, args;kw...) = Term{symtype(x)}(f, args)

    rules = [@rule (~a::isbp) ^ (~b::(x-> x isa Integer)) =>
             BoundedDegreeMonomial(((~a).p)^(~b), (~a).coeff ^ (~b), ~b > deg)
             @rule *(~~x) => mul_bounded(~~x, deg)]


    expr′ = Postwalk(Chain(rules), similarterm=st)(expr′)

    # Step 3: every term now has at most one valid monomial -- find the coefficient of it.
    expr′ = Postwalk(Chain([@rule *(~~x, ~a::isbp, ~~y) =>
                                BoundedDegreeMonomial((~a).p,
                                                      (~a).coeff * _mul(~~x) * _mul(~~y),
                                                      (~a).overdegree)]),
                     similarterm=st)(expr′)

    # Step 4: Bifurcate polynomial and nonlinear parts:

    #=
    if isbp(expr′)
        return Dict(expr′.p => expr.coeff), 0
    elseif istree(expr) && operation(expr) == (+)
        args = unsorted_arguments(expr)
        nls = filter(!isbp, args)
        nl = clean_sum(nls)
        filter(isbp, args), nl
    end
    =#
end

_mul(x) = isempty(x) ? 1 : isone(length(x)) ? first(x) : prod(x)

function mul_bounded(xs, deg)
    length(xs) == 1 ?
        first(xs) :
        reduce((x,y) -> mul(x,y, deg), xs)
end

pdegree(x::Mul) = sum(values(x.dict))
pdegree(x::Sym) = 1
pdegree(x::Pow) = x.exp
pdegree(x::Number) = 0

function mul(a::BoundedDegreeMonomial, b::BoundedDegreeMonomial, deg)
    if a.overdegree || b.overdegree || pdegree(a.p) + pdegree(b.p) > deg
        BoundedDegreeMonomial(1, a.p * b.p * a.coeff * b.coeff, true)
    else
        BoundedDegreeMonomial(a.p * b.p, a.coeff * b.coeff, false)
    end
end

mul(a::BoundedDegreeMonomial, b, deg) = BoundedDegreeMonomial(a.p, a.coeff * b, false)
mul(a, b::BoundedDegreeMonomial, deg) = BoundedDegreeMonomial(b.p, a * b.coeff, false)

function mark_vars(expr, vars)
    if !istree(expr)
        if any(isequal(expr), vars)
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
