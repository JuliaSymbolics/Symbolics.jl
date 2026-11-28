SymbolicUtils.promote_symtype(::typeof(imag), ::Type{Complex{T}}) where {T} = T
Base.promote_rule(::Type{Complex{T}}, ::Type{S}) where {T<:Real, S<:Num} =  Complex{S} # 283
Base.promote_rule(::Type{Complex{T}}, ::Type{Num}) where {T <: Real} = Complex{Num}

is_wrapper_type(::Type{Complex{Num}}) = true
has_symwrapper(::Type{<:Complex{T}}) where {T<:Real} = true
wraps_type(::Type{Complex{Num}}) = Complex{Real}
iswrapped(::Complex{Num}) = true
function wrapper_type(::Type{Complex{T}}) where T
    Symbolics.has_symwrapper(T) ? Complex{wrapper_type(T)} : Complex{T}
end

const BuiltinComplexInteger = Union{
    Bool, Int8, Int16, Int32, Int64, Int128, UInt8, UInt16, UInt32, UInt64, UInt128, BigInt,
}
const BuiltinComplexComponent = Union{
    BuiltinComplexInteger, Rational{<:BuiltinComplexInteger}, Float16, Float32, Float64, BigFloat,
}

# Cross-type equality for a custom Real may promote both values and hide rounding.
exact_complex(r::Real, i::Real) = Complex{typejoin(typeof(r), typeof(i))}(r, i)

function exact_complex(r::BuiltinComplexComponent, i::BuiltinComplexComponent)
    T = promote_type(typeof(r), typeof(i))
    z = try
        Complex{T}(r, i)
    catch err
        err isa InexactError || rethrow()
        return Complex{Real}(r, i)
    end
    return isequal(real(z), r) && isequal(imag(z), i) ? z : Complex{Real}(r, i)
end

function SymbolicUtils.unwrap(a::Complex{<:Num})
    re, img = real(unwrap(real(a))), real(unwrap(imag(a)))
    if SymbolicUtils.isconst(re) && SymbolicUtils.isconst(img)
        return Const{VartypeT}(exact_complex(unwrap_const(re), unwrap_const(img)))
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

const IM = Sym{VartypeT}(:im; type = Number)

function Base.show(io::IO, a::Complex{Num})
    rr = unwrap(real(a))
    ii = unwrap(imag(a))

    if iscall(rr) && (operation(rr) === real) &&
        iscall(ii) && (operation(ii) === imag) &&
        isequal(arguments(rr)[1], arguments(ii)[1])

        return print(io, arguments(rr)[1])
    end

    show(io, real(a) + IM * imag(a))
end

function (s::SymbolicUtils.Substituter)(x::Complex{Num})
    Complex{Num}(s(real(x)), s(imag(x)))
end
