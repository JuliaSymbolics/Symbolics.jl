@metatheory_init ()

using SymbolicUtils.Rewriters

TermInterface.isterm(t::Type{<:Sym}) = false
TermInterface.isterm(t::Type{<:Symbolic}) = true

TermInterface.gethead(t::Symbolic) = :call 
TermInterface.gethead(t::Sym) = t
TermInterface.getargs(t::Symbolic) = [operation(t), arguments(t)...]
TermInterface.arity(t::Symbolic) = length(arguments(t))

function TermInterface.similarterm(x::Type{<:Symbolic{T}}, head, args; metadata=nothing) where T
    @assert head == :call
    Term{T}(args[1], args[2:end])
end

"""
Converts `Mul` and `Add` to `Term`.
"""
function toterm(t::Mul{T}) where T
    args = []
    push!(args, t.coeff)
    for (k, deg) in t.dict
        push!(args, deg == 1 ? k : Term{T}(^, [k, deg]))
    end
    Term{T}(*, args)
end

function toterm(t::Add{T}) where T
    args = []
    for (k, coeff) in t.dict
        push!(args, coeff == 1 ? k : Term{T}(*, [coeff, k]))
    end
    Term{T}(+, args)
end

function toterm(t::Pow{T}) where T
    Term{T}(^, [t.base, t.exp])
end

toterm(t) = t

"""
Binarizes `Term`s with n-ary operations
"""
function unflatten(t::Symbolic{T}) where{T}
    # TODO change to isterm after PR
    if SymbolicUtils.istree(t)
        f = operation(t)
        if f == (+) || f == (*) || f == (-)  # check out for other binary ops TODO
            a = arguments(t)
            return foldl((x,y) -> Term{T}(f, [x, y]), a)
        end
    end
    return t
end

unflatten(t) = t

function preprocess(t)
    Chain([Postwalk(toterm), Postwalk(unflatten)])(t) 
end

# TODO turn back a tree of `Term` into `Mul`, `Add` and `Pow` types
#function rebuild(x::Term{T})
#    if operation(T) == * 
#        nothing
#    end
#end 

rebuild(x) = x

"""
Equational rewrite rules for optimizing expressions
"""
opt_theory = @methodtheory begin
    a + b == b + a
    a * b == b * a
    a * x + a * y == a*(x+y)
    -1 * a == -a
    a + (-1 * b) == a - b
    x^-1 == 1/x 
    1/x * a == a/x
    # fraction rules 
    # (a/b) + (c/b) => (a+c)*(1/b)
    # trig functions
    sin(x)/cos(x) == tan(x)
    cos(x)/sin(x) == cot(x)
    sin(x)^2 + cos(x)^2 => 1
    sin(2a) == 2sin(a)cos(a)
end

"""
Approximation of costs of operators in number 
of CPU cycles required for the numerical computation

See 
 * https://latkin.org/blog/2014/11/09/a-simple-benchmark-of-various-math-operations/
 * https://streamhpc.com/blog/2012-07-16/how-expensive-is-an-operation-on-a-cpu/
 * https://github.com/triscale-innov/GFlops.jl
"""
const op_costs = Dict(
    (+)     => 1,
    (-)     => 1,
    abs     => 2,
    (*)     => 3,
    exp     => 18,
    (/)     => 24,
    (^)     => 100,
    log1p   => 124,
    deg2rad => 125,
    rad2deg => 125,
    acos    => 127,
    asind   => 128,
    acsch   => 133,
    sin     => 134,
    cos     => 134,
    atan    => 135,
    tan     => 156,
)
# TODO some operator costs are in FLOP and not in cycles!!

function costfun(n::ENode, g::EGraph, an)
    arity(n) == 0 && return get(op_costs, n.head, 1)

    if !(n.head == :call)
        return 1000000000
    end
    cost = 0

    for id ∈ n.args
        eclass = geteclass(g, id)
        !hasdata(eclass, an) && (cost += Inf; break)
        cost += last(getdata(eclass, an))
    end
    cost
end

function optimize(ex; params=SaturationParams(timeout=20))
    prex = preprocess(unwrap(ex))
    g = EGraph()
    settermtype!(g, Term{symtype(ex), Any})
    ec, _ = addexpr!(g, prex)
    g.root = ec.id
    display(g.classes); println()
    saturate!(g, opt_theory, params)
    extract!(g, costfun) # --> "term" head args
end
