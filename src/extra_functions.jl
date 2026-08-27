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
    T1 == Num && T2 == AbstractArray && continue
    @eval function Base.in(x::$T1, y::$T2)
        return in(unwrap(x), unwrap(y))
    end
end
Base.in(x::Num, y::AbstractArray) = SymbolicUtils.term(in, unwrap(x), y)
Base.in(x::Num, y::Array) = SymbolicUtils.term(in, unwrap(x), y)
Base.in(x::Num, y::SparseArrays.AbstractSparseArray) = SymbolicUtils.term(in, unwrap(x), y)
Base.in(x::Num, y::AbstractRange{Num}) = SymbolicUtils.term(in, unwrap(x), y)
Base.in(x::Num, y::AbstractRange{<:Real}) = SymbolicUtils.term(in, unwrap(x), y)
Base.in(x::Num, y::AbstractRange{<:Integer}) = SymbolicUtils.term(in, unwrap(x), y)
function Base.in(x::Num, y::AbstractRange{T}) where {T >: Num}
    return SymbolicUtils.term(in, unwrap(x), y)
end
Base.in(x::Num, y::_StaticArray) = SymbolicUtils.term(in, unwrap(x), y)
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
Base.intersect(a::AbstractRange, b::Arr{<:Any, 1}) = intersect(a, unwrap(b))
Base.intersect(a::Arr{<:Any, 1}, b::AbstractRange) = intersect(unwrap(a), b)

LinearAlgebra.norm(x::Num, p::Real) = abs(x)

for f in (
        SpecialFunctions.besseli, SpecialFunctions.besselj,
        SpecialFunctions.besselk, SpecialFunctions.bessely,
    )
    @eval begin
        (f::$(typeof(f)))(a::Num, b::AbstractFloat) = invoke(f, Tuple{Num, Real}, a, b)
        (f::$(typeof(f)))(a::Num, b::Complex) = invoke(f, Tuple{Num, Number}, a, b)
    end
end
for T in (Signed, Float32, Float64)
    @eval Base.copysign(a::$T, b::Num) = invoke(copysign, Tuple{Real, Num}, a, b)
end
SpecialFunctions.polygamma(a::Integer, b::Num) =
    invoke(SpecialFunctions.polygamma, Tuple{Real, Num}, a, b)

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
