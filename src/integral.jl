# basic integral struct with upper bound and lower bound.
struct Integral<: Function
    x
    ub
    lb
    Integral(x,ub,lb) = new(Symbolics.value(x), Symbolics.value(ub), Symbolics.value(lb))
end

(I::Integral)(x) = Term{SymbolicUtils.symtype(x)}(I, [x])
(I::Integral)(x::Num) = Num(I(Symbolics.value(x)))
SymbolicUtils.promote_symtype(::Integral, x) = x

function Base.show(io::IO, I::Integral)
    print(io, "Integral(", I.x, ")")
    print(io,"upper_bound: ")
    show(io,I.ub)
    println(io , "lower_bound")
    show(io , I.lb)
end

Base.:(==)(I1::Integral, I2::Integral) = isequal(I1.x, I2.x)
