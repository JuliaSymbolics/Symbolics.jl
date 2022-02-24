using Symbolics
using ChainRules
using ChainRulesCore
using Symbolics: unwrap, wrap
using SymbolicUtils
using SymbolicUtils: unsorted_arguments, issym
using SymbolicUtils.Rewriters

using SparseArrays
struct Dual <: Real
    primal::Real
    partial::AbstractVector # View(SparseMarixCSC, ...)
end

function diff(f, x)
    duals = map(enumerate(xs)) do (i, x,)
        ps = spzeros(n)
        ps[i] = 1
        Dual(x, ps)
    end

    y = unwrap(substitute(f, Dict(x => Dual(x,1)), fold=true))
    y.partial
end

function jac(f, xs)
    n = length(xs)
    duals = map(enumerate(xs)) do (i, x,)
        ps = spzeros(n)
        ps[i] = 1
        Dual(x, ps)
    end

    ys = unwrap.(substitute.(f, (Dict(xs .=> duals),), fold=true))
    hcat(map(y->y.partial, ys)...)'
end

function df(f, xs...)
    ChainRules.frule((ChainRules.NoTangent(),
    map(x->(x isa Dual ? x.partial : 0), xs)...),
    f, map(x->x isa Dual ? x.primal : x, xs)...)
end

SymbolicUtils.@number_methods Dual Dual(df(f, a)...) Dual(df(f, a, b)...)

Base.:(^)(d::Dual, i::Integer) = Dual(df(^, d, i)...)

struct DualArray{T,N} <: AbstractArray{T,N}
    primal::AbstractArray{T,N}
    partial::AbstractArray # View(SparseMarixCSC, ...)
end

Base.size(d::DualArray) = size(d.primal)
Dual(x::AbstractArray, dx) = DualArray(x, dx)
Dual(x::AbstractArray) = DualArray(x, onehot(length(x)))

function gradient(fx)
    fx, jvp = rrule(operation(fx), arguments(fx)...)
end

function onehot(n)
    sparse(1:n,1:n,ones(Num, n))
end
