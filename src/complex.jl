include("symbolic_number.jl")

const SymbolicScalar = Union{Num, SymbolicNumber, Complex{Num}}

SymbolicUtils.promote_symtype(::typeof(imag), ::Type{Complex{T}}) where {T} = T
Base.promote_rule(::Type{Complex{T}}, ::Type{Num}) where {T <: Real} = SymbolicNumber
Base.promote_rule(::Type{Num}, ::Type{Complex{T}}) where {T <: Real} = SymbolicNumber

# Wrapper arithmetic narrows at runtime, so inference cannot pick a concrete eltype.
Base.promote_op(
    ::typeof(LinearAlgebra.matprod), ::Type{SymbolicNumber}, ::Type{SymbolicNumber}
) = SymbolicNumber
Base.promote_op(
    ::typeof(LinearAlgebra.matprod), ::Type{SymbolicNumber}, ::Type{T}
) where {T <: Number} = SymbolicNumber
Base.promote_op(
    ::typeof(LinearAlgebra.matprod), ::Type{T}, ::Type{SymbolicNumber}
) where {T <: Number} = SymbolicNumber

for C in (Complex, Complex{Bool})
    @eval begin
        Base.:+(x::Num, z::$C) = wrap(unwrap(x) + z)
        Base.:+(z::$C, x::Num) = wrap(z + unwrap(x))
        Base.:-(x::Num, z::$C) = wrap(unwrap(x) - z)
        Base.:-(z::$C, x::Num) = wrap(z - unwrap(x))
        Base.:*(x::Num, z::$C) = wrap(unwrap(x) * z)
        Base.:*(z::$C, x::Num) = wrap(z * unwrap(x))
        Base.:/(x::Num, z::$C) = wrap(unwrap(x) / z)
        Base.:/(z::$C, x::Num) = wrap(z / unwrap(x))
    end
end

const ConcreteReal = Union{AbstractFloat, Integer, Rational, AbstractIrrational}
Base.:^(a::Complex{T}, b::Num) where {T <: ConcreteReal} = wrap(term(^, a, unwrap(b)))

is_wrapper_type(::Type{Complex{Num}}) = true
wraps_type(::Type{Complex{Num}}) = Complex{Real}
iswrapped(::Complex{Num}) = true

function SymbolicUtils.unwrap(a::Complex{<:Num})
    re, img = unwrap(real(a)), unwrap(imag(a))
    if SymbolicUtils.isconst(re) && SymbolicUtils.isconst(img)
        return Const{VartypeT}(complex(unwrap_const(re), unwrap_const(img)))
    end
    if iscall(re) && operation(re) === real && iscall(img) && operation(img) === imag && isequal(arguments(re)[1], arguments(img)[1])
        return arguments(re)[1]
    end
    sT = promote_type(symtype(re), symtype(img))
    return Term{VartypeT}(complex, SymbolicUtils.ArgsT{vartype(re)}((re, img)); type = Complex{sT}, shape = SymbolicUtils.ShapeVecT())
end

SymbolicUtils.infer_vartype(::Type{Complex{Num}}) = VartypeT

function Base.Complex{Num}(x::BasicSymbolic{VartypeT})
    Complex{Num}(wrap(real(x)), wrap(imag(x)))
end

function Base.show(io::IO, a::Complex{Num})
    rr = unwrap(real(a))
    ii = unwrap(imag(a))

    if iscall(rr) && (operation(rr) === real) &&
        iscall(ii) && (operation(ii) === imag) &&
        isequal(arguments(rr)[1], arguments(ii)[1])

        return print(io, arguments(rr)[1])
    end

    return show(io, real(a) + im * imag(a))
end

function (s::SymbolicUtils.Substituter)(x::Complex{Num})
    Complex{Num}(s(real(x)), s(imag(x)))
end

# Base's `Complex` methods branch on component values, which symbolic parts cannot do.
function Base.exp(z::Complex{Num})
    a, b = reim(z)
    ea = exp(a)
    return Complex(ea * cos(b), ea * sin(b))
end
function Base.sin(z::Complex{Num})
    a, b = reim(z)
    return Complex(sin(a) * cosh(b), cos(a) * sinh(b))
end
function Base.cos(z::Complex{Num})
    a, b = reim(z)
    return Complex(cos(a) * cosh(b), -sin(a) * sinh(b))
end
function Base.log(z::Complex{Num})
    a, b = reim(z)
    r = sqrt(a^2 + b^2)
    return Complex(log(r), atan(b, a))
end
function Base.sqrt(z::Complex{Num})
    a, b = reim(z)
    r = sqrt(a^2 + b^2)
    θ = atan(b, a) / 2
    sr = sqrt(r)
    return Complex(sr * cos(θ), sr * sin(θ))
end

# Principal branch.
function _cartesian_pow(z::Complex{Num}, p::Real)
    a, b = reim(z)
    r = sqrt(a^2 + b^2)
    θ = atan(b, a)
    rp = r^p
    pθ = p * θ
    return Complex(rp * cos(pθ), rp * sin(pθ))
end

Base.:^(z::Complex{Num}, p::AbstractFloat) = _cartesian_pow(z, p)
Base.:^(z::Complex{Num}, p::Rational) = _cartesian_pow(z, p)
Base.:^(z::Complex{Num}, p::Num) = _cartesian_pow(z, p)
