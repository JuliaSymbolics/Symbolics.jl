SymbolicUtils.promote_symtype(::typeof(imag), ::Type{Complex{T}}) where {T} = T
Base.promote_rule(::Type{Complex{T}}, ::Type{S}) where {T<:Real, S<:Num} =  Complex{S} # 283

is_wrapper_type(::Type{Complex{Num}}) = true
has_symwrapper(::Type{<:Complex{T}}) where {T<:Real} = true
wraps_type(::Type{Complex{Num}}) = Complex{Real}
iswrapped(::Complex{Num}) = true
function wrapper_type(::Type{Complex{T}}) where T
    Symbolics.has_symwrapper(T) ? Complex{wrapper_type(T)} : Complex{T}
end

function SymbolicUtils.unwrap(a::Complex{<:Num})
    re, img = unwrap(real(a)), unwrap(imag(a))
    sT = promote_type(symtype(re), symtype(img))
    return Term{vartype(re)}(complex, SymbolicUtils.ArgsT{vartype(re)}((re, img)); type = Complex{sT}, shape = SymbolicUtils.ShapeVecT())
end
