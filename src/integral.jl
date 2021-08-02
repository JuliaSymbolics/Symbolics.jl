# basic integral struct with upper bound and lower bound.
struct Integral{X, T <: Symbolics.VarDomainPairing} <: Function
    x
    domain::T
    Integral(x,domain) = new{typeof(x),typeof(domain)}(Symbolics.value.(x), domain)
end

(I::Integral)(x) = Term{SymbolicUtils.symtype(x)}(I, [x])
(I::Integral)(x::Num) = Num(I(Symbolics.value(x)))
SymbolicUtils.promote_symtype(::Integral, x) = x

function Base.show(io::IO, I::Integral)
    print(io, "Integral(", I.x, ", ", I.domain, ")")
end

Base.:(==)(I1::Integral, I2::Integral) = convert(Bool, simplify(isequal(I1.x, I2.x) && isequal(I1.domain, I2.domain)))
