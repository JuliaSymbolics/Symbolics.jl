include("symbolic_number.jl")

const SymbolicScalar = Union{Num, SymbolicNumber, Complex{Num}}

SymbolicUtils.promote_symtype(::typeof(imag), ::Type{Complex{T}}) where {T} = T
Base.promote_rule(::Type{Complex{T}}, ::Type{Num}) where {T <: Real} = SymbolicNumber
Base.promote_rule(::Type{Num}, ::Type{Complex{T}}) where {T <: Real} = SymbolicNumber

function _symbolic_matmul(A, B::AbstractMatrix)
    return mul!(Matrix{SymbolicNumber}(undef, size(A, 1), size(B, 2)), A, B)
end

function _symbolic_matmul(A, b::AbstractVector)
    return mul!(Vector{SymbolicNumber}(undef, size(A, 1)), A, b)
end

Base.:*(A::StridedMatrix{SymbolicNumber}, B::StridedVecOrMat{SymbolicNumber}) = _symbolic_matmul(A, B)
Base.:*(A::StridedMatrix{SymbolicNumber}, B::StridedVecOrMat{<:Number}) = _symbolic_matmul(A, B)
Base.:*(A::StridedMatrix{<:Number}, B::StridedVecOrMat{SymbolicNumber}) = _symbolic_matmul(A, B)

Base.:*(A::SparseArrays.SparseMatrixCSC{SymbolicNumber}, B::StridedVecOrMat{SymbolicNumber}) = _symbolic_matmul(A, B)
Base.:*(A::SparseArrays.SparseMatrixCSC{SymbolicNumber}, B::StridedVecOrMat{<:Number}) = _symbolic_matmul(A, B)
Base.:*(A::SparseArrays.SparseMatrixCSC{<:Number}, B::StridedVecOrMat{SymbolicNumber}) = _symbolic_matmul(A, B)

function _cartesian_parts(x)
    symtype(x) <: Real && return (x, 0)
    if iscall(x)
        op = operation(x)
        if op === (+)
            re, im_part = 0, 0
            for arg in arguments(x)
                ar, ai = _cartesian_parts(arg)
                re += ar
                im_part += ai
            end
            return re, im_part
        elseif op === (*)
            re, im_part = 1, 0
            for arg in arguments(x)
                ar, ai = _cartesian_parts(arg)
                re, im_part = re * ar - im_part * ai, re * ai + im_part * ar
            end
            return re, im_part
        end
    end
    return real(x), imag(x)
end

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
    return Complex{Num}(wrap(real(x)), wrap(imag(x)))
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
    return Complex{Num}(s(real(x)), s(imag(x)))
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
function _cartesian_call(f, z::Complex{Num}, args...)
    raw = unwrap(z)
    if SymbolicUtils.isconst(raw) && all(SymbolicUtils.isconst, args)
        result = f(unwrap_const(raw), map(unwrap_const, args)...)
        return Complex{Num}(result)
    end
    return Complex{Num}(term(f, raw, args...; type = Complex{Real}))
end

Base.log(z::Complex{Num}) = _cartesian_call(log, z)
Base.sqrt(z::Complex{Num}) = _cartesian_call(sqrt, z)
Base.:^(z::Complex{Num}, p::AbstractFloat) = _cartesian_call(^, z, unwrap(p))
Base.:^(z::Complex{Num}, p::Rational) = _cartesian_call(^, z, unwrap(p))
Base.:^(z::Complex{Num}, p::Num) = _cartesian_call(^, z, unwrap(p))
