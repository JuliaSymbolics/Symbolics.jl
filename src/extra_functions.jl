for (T1, T2) in Iterators.product(Iterators.repeated([Number, BasicSymbolic{VartypeT}, Num], 2)...)
    if T1 != Num && T2 != Num
        continue
    end
    @eval function Base.binomial(a::$T1, b::$T2)
        binomial(unwrap(a), unwrap(b))
    end
end

derivative(::typeof(sign), args::NTuple{1,Any}, ::Val{1}) = 0

derivative(::typeof(signbit), args::NTuple{1,Any}, ::Val{1}) = 0
derivative(::typeof(abs), args::NTuple{1,Any}, ::Val{1}) = ifelse(signbit(args[1]),-one(args[1]),one(args[1]))

function derivative(::typeof(min), args::NTuple{2,Any}, ::Val{1})
    x, y = args
    ifelse(x < y, one(x), zero(x))
end
function derivative(::typeof(min), args::NTuple{2,Any}, ::Val{2})
    x, y = args
    ifelse(x < y, zero(y), one(y))
end
function derivative(::typeof(max), args::NTuple{2,Any}, ::Val{1})
    x, y = args
    ifelse(x > y, one(x), zero(x))
end
function derivative(::typeof(max), args::NTuple{2,Any}, ::Val{2})
    x, y = args
    ifelse(x > y, zero(y), one(y))
end

function derivative(::Union{typeof(ceil),typeof(floor),typeof(factorial)}, args::NTuple{1,Any}, ::Val{1})
    zero(args[1])
end

@register_symbolic Base.rand(x)
@register_symbolic Base.randn(x)

for (T1, T2, T3) in Iterators.product(Iterators.repeated((Num, BasicSymbolic{VartypeT}, Number), 3)...)
    if T1 != Num && T2 != Num && T3 != Num
        continue
    end
    @eval function Base.clamp(a::$T1, b::$T2, c::$T3)
        wrap(clamp(unwrap(a), unwrap(b), unwrap(c)))
    end
end

function derivative(::typeof(Base.clamp), args::NTuple{3, Any}, ::Val{1})
    x, l, h = args
    T = promote_type(symtype(x), symtype(l), symtype(h))
    z = zero(T)
    o = one(T)
    ifelse(x<l, z, ifelse(x>h, z, o))
end

for T1 in [Real, Num, BasicSymbolic{VartypeT}], T2 in [AbstractArray, Arr, BasicSymbolic{VartypeT}]
    if T1 != Num && T2 != Arr
        continue
    end
    @eval function Base.in(x::$T1, y::$T2)
        return in(unwrap(x), unwrap(y))
    end
end
for (T1, T2) in Iterators.product(Iterators.repeated([AbstractArray, Arr, BasicSymbolic{VartypeT}], 2)...)
    if T1 != Arr && T2 != Arr
        continue
    end
    @eval function Base.union(a::$T1, b::$T2)
        union(unwrap(a), unwrap(b))
    end
    @eval function Base.intersect(a::$T1, b::$T2)
        intersect(unwrap(a), unwrap(b))
    end
    @eval function Base.issubset(a::$T1, b::$T2)
        issubset(unwrap(a), unwrap(b))
    end
end

LinearAlgebra.norm(x::Num, p::Real) = abs(x)

derivative(::typeof(<), ::NTuple{2, Any}, ::Val{i}) where {i} = 0
derivative(::typeof(<=), ::NTuple{2, Any}, ::Val{i}) where {i} = 0
derivative(::typeof(>), ::NTuple{2, Any}, ::Val{i}) where {i} = 0
derivative(::typeof(>=), ::NTuple{2, Any}, ::Val{i}) where {i} = 0
derivative(::typeof(==), ::NTuple{2, Any}, ::Val{i}) where {i} = 0
derivative(::typeof(!=), ::NTuple{2, Any}, ::Val{i}) where {i} = 0

derivative(::typeof(expinti), args::NTuple{1,Any}, ::Val{1}) = exp(args[1])/args[1]
