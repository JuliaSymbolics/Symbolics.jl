# basic integral struct with upper bound and lower bound.
"""
    Integral(domain)

Defines an Integral operator I(ex) which represents the integral of `I` of the
expression `ex` over the `domain`. Note that the `domain` must be a 
`Symbolics.VarDomainPairing` where the chosen variable is the variable being
integrated over, i.e. `Integral(x in domain)` means that `I` is the integral
operator with respect to `dx`.

## Examples

```julia
I1 = Integral(x in ClosedInterval(1, 5))
I2 = Integral(x in ClosedInterval(1, 5))

@variables a b
I = Integral(x in ClosedInterval(a, b))
@test isequal(I(0), 0)
@test isequal(I(2), 2*(b -a))
```
"""
struct Integral
    domain::VarDomainPairing
end

Base.nameof(::Integral) = :Integral

const AnyInterval{T} = DomainSets.Interval{A, B, T} where {A, B}

function (I::Integral)(x::Union{Rational, AbstractIrrational, AbstractFloat, Integer})
    domain = I.domain.domain
    # Despite this not being a concrete type, `endpoints` only needs to access fields
    # of `Interval` so this ends up being type-stable. There doesn't seem to be a way
    # to test this and ensure it remains type-stable.
    if domain isa AnyInterval{Num}
        a, b = unwrap.(DomainSets.endpoints(domain))
        return Num((b - a) * x)
    elseif domain isa AnyInterval{SymbolicT}
        a, b = DomainSets.endpoints(domain)
        return Num((b - a) * x)
    elseif domain isa AnyInterval{Int}
        a, b = DomainSets.endpoints(domain)
        return Num((b - a) * x)
    elseif domain isa AnyInterval{Float64}
        a, b = DomainSets.endpoints(domain)
        return Num((b - a) * x)
    else
        # ::NTuple{2, Any} avoids `indexed_iterate` dynamic dispatching
        a, b = DomainSets.endpoints(domain)::NTuple{2, Any}
        # SConst avoids `*` dynamic dispatching
        return Num(SConst(b - a) * x)
    end
end
(I::Integral)(x::Complex) = Complex{Num}(Num(I(unwrap(real(x)))), Num(I(unwrap(imag(x)))))
function (I::Integral)(x)
    return Term{VartypeT}(
        I, SArgsT((x,));
        type = SymbolicUtils.symtype(x), shape = SymbolicUtils.shape(x)
    )
end
(I::Integral)(x::Num) = Num(I(unwrap(x)))
SymbolicUtils.promote_symtype(::Integral, T::SymbolicUtils.TypeT) = T
SymbolicUtils.promote_shape(::Integral, @nospecialize(sh::SymbolicUtils.ShapeT)) = sh

function Base.show(io::IO, I::Integral)
    print(io, "Integral(", I.domain.variables, ", ", I.domain.domain, ")")
end

Base.:(==)(I1::Integral, I2::Integral) = convert(Bool, simplify(isequal(I1.domain, I2.domain)))
