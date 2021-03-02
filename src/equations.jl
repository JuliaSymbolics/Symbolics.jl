"""
$(TYPEDEF)

An equality relationship between two expressions.

# Fields
$(FIELDS)
"""
struct Equation
    """The expression on the left-hand side of the equation."""
    lhs
    """The expression on the right-hand side of the equation."""
    rhs
end
Base.:(==)(a::Equation, b::Equation) = all(isequal.((a.lhs, a.rhs), (b.lhs, b.rhs)))
Base.hash(a::Equation, salt::UInt) = hash(a.lhs, hash(a.rhs, salt))

Base.show(io::IO, eq::Equation) = print(io, eq.lhs, " ~ ", eq.rhs)

SymbolicUtils.simplify(x::Equation; kw...) = simplify(x.lhs; kw...) ~ simplify(x.rhs; kw...)

lhss(xs) = map(x->x.lhs, xs)
rhss(xs) = map(x->x.rhs, xs)

"""
$(TYPEDSIGNATURES)

Create an [`Equation`](@ref) out of two [`Num`](@ref) instances, or an
`Num` and a `Number`.

# Examples

```jldoctest
julia> using Symbolics

julia> @variables x y;

julia> x ~ y
x ~ y

julia> x - y ~ 0
x - y ~ 0
```
"""
Base.:~(lhs::Num, rhs::Num) = Equation(value(lhs), value(rhs))
Base.:~(lhs::Num, rhs::Number    ) = Equation(value(lhs), value(rhs))
Base.:~(lhs::Number    , rhs::Num) = Equation(value(lhs), value(rhs))
Base.:~(lhs::Symbolic, rhs::Symbolic) = Equation(value(lhs), value(rhs))
Base.:~(lhs::Symbolic, rhs::Any    ) = Equation(value(lhs), value(rhs))
Base.:~(lhs::Any, rhs::Symbolic    ) = Equation(value(lhs), value(rhs))
for T in [:Num, :Complex, :Number], S in [:Num, :Complex, :Number]
    (T != :Complex && S != :Complex) && continue
    @eval Base.:~(a::$T, b::$S) = [
      (isa(value(real(a)), Number) && isa(value(real(b)), Number)) ? [] : value(real(a)) ~ value(real(b));
      (isa(value(imag(a)), Number) && isa(value(imag(b)), Number)) ? [] : value(imag(a)) ~ value(imag(b))
    ]
end

struct ConstrainedEquation
  constraints
  eq
end

function expand_derivatives(eq::Equation, simplify=false)
    return Equation(expand_derivatives(eq.lhs, simplify), expand_derivatives(eq.rhs, simplify))
end
