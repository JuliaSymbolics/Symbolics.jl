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
function SymbolicUtils.substitute(x::Equation, rules; kw...)
    sub = substituter(rules)
    sub(x.lhs; kw...) ~ sub(x.rhs; kw...)
end

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

julia> @variables A[1:3, 1:3] B[1:3, 1:3];

julia> x ~ y
x ~ y

julia> x - y ~ 0
x - y ~ 0

julia> A .~ B
3×3 Array{Equation,2}:
 A₁ˏ₁ ~ B₁ˏ₁  A₁ˏ₂ ~ B₁ˏ₂  A₁ˏ₃ ~ B₁ˏ₃
 A₂ˏ₁ ~ B₂ˏ₁  A₂ˏ₂ ~ B₂ˏ₂  A₂ˏ₃ ~ B₂ˏ₃
 A₃ˏ₁ ~ B₃ˏ₁  A₃ˏ₂ ~ B₃ˏ₂  A₃ˏ₃ ~ B₃ˏ₃

julia> A .~ 3x
3×3 Array{Equation,2}:
 A₁ˏ₁ ~ 3x  A₁ˏ₂ ~ 3x  A₁ˏ₃ ~ 3x
 A₂ˏ₁ ~ 3x  A₂ˏ₂ ~ 3x  A₂ˏ₃ ~ 3x
 A₃ˏ₁ ~ 3x  A₃ˏ₂ ~ 3x  A₃ˏ₃ ~ 3x
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
    @eval Base.:~(a::$T, b::$S) = let ar = value(real(a)), br = value(real(b)),
                                      ai = value(imag(a)), bi = value(imag(b))
        if ar isa Number && br isa Number && ai isa number && bi isa Number
            error("Equation $a ~ $b does not contain any symbols")
        elseif ar isa Number && br isa Number
            ai ~ bi
        elseif ai isa Number && bi isa Number
            ar ~ br
        else
            [ar ~ br
            ai ~ bi]
        end
    end
end

struct ConstrainedEquation
  constraints
  eq
end

function expand_derivatives(eq::Equation, simplify=false)
    return Equation(expand_derivatives(eq.lhs, simplify), expand_derivatives(eq.rhs, simplify))
end
