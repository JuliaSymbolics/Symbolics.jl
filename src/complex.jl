SymbolicUtils.promote_symtype(::typeof(imag), ::Type{Complex{T}}) where {T} = T
Base.promote_rule(::Type{Complex{T}}, ::Type{S}) where {T <: Real, S <: Num} = Complex{S} # 283

is_wrapper_type(::Type{Complex{Num}}) = true
has_symwrapper(::Type{<:Complex{T}}) where {T <: Real} = true
wraps_type(::Type{Complex{Num}}) = Complex{Real}
iswrapped(::Complex{Num}) = true
function wrapper_type(::Type{Complex{T}}) where {T}
    Symbolics.has_symwrapper(T) ? Complex{wrapper_type(T)} : Complex{T}
end

function SymbolicUtils.unwrap(a::Complex{<:Num})
    re, img = unwrap(real(a)), unwrap(imag(a))
    if SymbolicUtils.isconst(re) && SymbolicUtils.isconst(img)
        return Const{VartypeT}(complex(unwrap_const(re), unwrap_const(img)))
    end
    if iscall(re) && operation(re) === real && iscall(img) && operation(img) === imag &&
       isequal(arguments(re)[1], arguments(img)[1])
        return arguments(re)[1]
    end
    sT = promote_type(symtype(re), symtype(img))
    return Term{VartypeT}(complex, SymbolicUtils.ArgsT{vartype(re)}((re, img));
        type = Complex{sT}, shape = SymbolicUtils.ShapeVecT())
end

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
    re_sub = s(real(x))
    im_sub = s(imag(x))
    # Unwrap to check if the substituted values are complex constants
    re_unwrapped = unwrap(re_sub)
    im_unwrapped = unwrap(im_sub)
    # If both parts are constants, we can evaluate the full complex expression
    if SymbolicUtils.isconst(re_unwrapped) && SymbolicUtils.isconst(im_unwrapped)
        re_val = SymbolicUtils.unwrap_const(re_unwrapped)
        im_val = SymbolicUtils.unwrap_const(im_unwrapped)
        # Properly handle complex arithmetic: (a + b*im) where a, b can be complex
        # (a_re + a_im*im) + (b_re + b_im*im)*im = (a_re - b_im) + (a_im + b_re)*im
        result = re_val + im_val * im
        return Complex{Num}(wrap(Const{VartypeT}(real(result))), wrap(Const{VartypeT}(imag(result))))
    end
    Complex{Num}(re_sub, im_sub)
end
