for (T1, T2) in Iterators.product([Number, BasicSymbolic{VartypeT}, Num], [Integer, BasicSymbolic{VartypeT}, Num])
    if T1 != Num && T2 != Num
        continue
    end
    @eval function Base.binomial(a::$T1, b::$T2)
        binomial(unwrap(a), unwrap(b))
    end
end

@register_derivative sign(x) 1 COMMON_ZERO
@register_derivative signbit(x) 1 COMMON_ZERO
@register_derivative abs(x) 1 ifelse(signbit(x),-one(x),one(x))
@register_derivative min(x, y) 1 ifelse(x < y, one(x), zero(x))
@register_derivative min(x, y) 2 ifelse(x < y, zero(y), one(y))
@register_derivative max(x, y) 1 ifelse(x > y, one(x), zero(x))
@register_derivative max(x, y) 2 ifelse(x > y, zero(y), one(y))
@register_derivative ceil(x) 1 COMMON_ZERO
@register_derivative floor(x) 1 COMMON_ZERO
@register_derivative factorial(x) 1 COMMON_ZERO

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

@register_derivative clamp(x, l, h) 1 ifelse(x < l, COMMON_ZERO, ifelse(x > h, COMMON_ZERO, COMMON_ONE))

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

# `Base.sincospi(::Real)` explicitly throws a `MethodError` so that subtypes of
# `Real` (like `Num`) have to opt in. Fall back to a tuple of `sinpi`/`cospi`,
# both of which are already registered.
Base.sincospi(x::Num) = (sinpi(x), cospi(x))

# `Base.sinpi(::Complex)` and `Base.cospi(::Complex)` branch on `isinteger` of
# the real/imaginary parts, so we need `isinteger(::Num)`. For a wrapped
# constant we can answer exactly; for a symbolic variable whose `symtype` is an
# `Integer` we know the value must be an integer; otherwise we conservatively
# return `false` (an unknown real is not assumed to be an integer).
function Base.isinteger(x::Num)
    val = unwrap(x)
    if SymbolicUtils.isconst(val)
        return isinteger(SymbolicUtils.unwrap_const(val))
    end
    return symtype(val) <: Integer
end

@register_derivative <(x, y) I COMMON_ZERO
@register_derivative <=(x, y) I COMMON_ZERO
@register_derivative >(x, y) I COMMON_ZERO
@register_derivative >=(x, y) I COMMON_ZERO
@register_derivative ==(x, y) I COMMON_ZERO
@register_derivative !=(x, y) I COMMON_ZERO
@register_derivative expinti(x) 1 exp(x) / x
