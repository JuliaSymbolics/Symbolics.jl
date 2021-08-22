# basic integral struct with upper bound and lower bound.
struct Integral{T <: Symbolics.VarDomainPairing} <: Function
    domain::T
    Integral(domain) = new{typeof(domain)}(domain)
end

(I::Integral)(x) = Term{SymbolicUtils.symtype(x)}(I, [x])
(I::Integral)(x::Num) = Num(I(Symbolics.value(x)))
SymbolicUtils.promote_symtype(::Integral, x) = x

function Base.show(io::IO, I::Integral)
    print(io, "Integral(", I.domain.variables, ", ", I.domain.domain, ")")
end

Base.:(==)(I1::Integral, I2::Integral) = convert(Bool, simplify(isequal(I1.domain, I2.domain)))
