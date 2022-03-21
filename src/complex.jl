abstract type AbstractComplexTerm{T} <: Symbolic{Complex{T}}
end

struct ComplexTerm{T} <: AbstractComplexTerm{T}
    re
    im
end

Base.imag(c::Symbolic{Complex{T}}) where {T} = term(imag, c)
SymbolicUtils.promote_type(::typeof(imag), ::Type{Complex{T}}) where {T} = T
Base.promote_rule(::Type{Complex{T}}, ::Type{S}) where {T<:Real, S<:Num} =  Complex{S} # 283

has_symwrapper(::Type{<:Complex{T}}) where {T<:Real} = true
wraps_type(::Type{Complex{Num}}) = Complex{Real}
iswrapped(::Complex{Num}) = true
function wrapper_type(::Type{Complex{T}}) where T
    Symbolics.has_symwrapper(T) ? Complex{wrapper_type(T)} : Complex{T}
end

TermInterface.symtype(a::ComplexTerm{T}) where T = Complex{T}
TermInterface.istree(a::ComplexTerm) = true
TermInterface.operation(a::ComplexTerm{T}) where T = Complex{T}
TermInterface.arguments(a::ComplexTerm) = [a.re, a.im]

function TermInterface.similarterm(t::ComplexTerm, f, args, symtype; metadata=nothing, exprhead=exprhead(t))
    if f <: Complex
        ComplexTerm{real(f)}(args...)
    else
        similarterm(first(args), f, args, symtype; metadata=metadata, exprhead=exprhead)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", a::ComplexTerm)
    print(io, "ComplexTerm(")
    show(io, mime, wrap(a))
    print(io, ")")
end

function Base.show(io::IO, a::Complex{Num})
    rr = unwrap(real(a))
    ii = unwrap(imag(a))

    if istree(rr) && (operation(rr) === real) &&
        istree(ii) && (operation(ii) === imag) &&
        isequal(arguments(rr)[1], arguments(ii)[1])

        return print(io, arguments(rr)[1])
    end

    i = Sym{Real}(:im)
    show(io, real(a) + i * imag(a))
end

function unwrap(a::Complex{<:Num})
    re, img = unwrap(real(a)), unwrap(imag(a))
    if re isa Real && img isa Real
        return re + im * img
    end
    T = promote_type(symtype(re), symtype(img))
    ComplexTerm{T}(re, img)
end
wrap(a::ComplexTerm) = Complex(wrap.(arguments(a))...)
wrap(a::Symbolic{<:Complex}) = Complex(wrap(real(a)), wrap(imag(a)))

SymbolicUtils.@number_methods(
                              ComplexTerm,
                              unwrap(f(wrap(a))),
                              unwrap(f(wrap(a), wrap(b))),
                             )

function Base.isequal(a::ComplexTerm{T}, b::ComplexTerm{S}) where {T,S}
    T === S && isequal(a.re, b.re) && isequal(a.im, b.im)
end

function Base.hash(a::ComplexTerm{T}, h::UInt) where T
    hash(hash(a.im, hash(a.re, hash(T, hash(h ⊻ 0x1af5d7582250ac4a)))))
end

Base.iszero(x::Complex{<:Num}) = iszero(real(x)) && iszero(imag(x))
Base.isone(x::Complex{<:Num}) = isone(real(x)) && iszero(imag(x))
_iszero(x::Complex{<:Num}) = _iszero(unwrap(x))
_isone(x::Complex{<:Num}) = _isone(unwrap(x))
