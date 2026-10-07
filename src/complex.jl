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

function SymbolicUtils.unwrap(a::Complex{<:Num})
    re, img = unwrap(real(a)), unwrap(imag(a))
    if SymbolicUtils.isconst(re) && SymbolicUtils.isconst(img)
        # `re`/`img` may wrap complex constants, for which `complex` would throw
        return Const{VartypeT}(unwrap_const(re) + im * unwrap_const(img))
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

# Split `x` into `(re, img)` with `x == re + im*img` and real-symtyped parts.
# e.g. `x => im*y` substituted into `1.7x` gives `1.7*complex(0, y)` with re/im `0, 1.7y`.
function _complex_reim(x::BasicSymbolic{VartypeT})
    if !(symtype(x) <: Real) && iscall(x)
        op = operation(x)
        if op === (+) || op === (*)
            args = arguments(x)
            re, img = _complex_reim(args[1])
            for i in 2:length(args)
                are, aim = _complex_reim(args[i])
                if op === (+)
                    re += are
                    img += aim
                else
                    re, img = re * are - img * aim, re * aim + img * are
                end
            end
            return re, img
        end
    end
    return real(x), imag(x)
end

function (s::SymbolicUtils.Substituter)(x::Complex{Num})
    val = get(SymbolicUtils.get_substitution_dict(s), x, nothing)
    val === nothing || return val
    are, aim = _complex_reim(unwrap(s(real(x))))
    re2, im2 = _complex_reim(unwrap(s(imag(x))))
    return Complex{Num}(wrap(are - im2), wrap(aim + re2))
end
