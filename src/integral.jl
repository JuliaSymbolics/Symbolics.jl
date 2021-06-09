# basic integral struct with upper bound and lower bound.
struct Integral{X, T <: Domain} <: Function
    x
    domain::T
    Integral(x,domain) = new{typeof(x),typeof(domain)}(Symbolics.value(x), domain)
end

(I::Integral)(x) = Term{SymbolicUtils.symtype(x)}(I, [x])
(I::Integral)(x::Num) = Num(I(Symbolics.value(x)))
SymbolicUtils.promote_symtype(::Integral, x) = x

function Base.show(io::IO, I::Integral)
    print(io, "Integral(", I.x, ", ", I.domain, ")")
end

Base.:(==)(I1::Integral, I2::Integral) = (isequal(I1.x, I2.x) && isequal(I1.domain, I2.domain))
(D::Differential)(I::Integral{X,T}) where{X,T} = I∘D

function replaceSym(a::Sym, b, O)
    if isa(O , Sym)
        if isequal(O , a)
            return b
        else
            return O
        end
    else
        args_replace = Vector{Symbolic{Real}}()
        args_ = arguments(O)
        for i in 1:length(args_)
            push!(args_replace, replaceSym(a , b , args_[i]))
        end
        return operation(O)(args_replace...)
    end
end
