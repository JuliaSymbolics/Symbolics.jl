const NAMESPACE_SEPARATOR = '₊'

hide_lhs(_) = false

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
    function Equation(lhs, rhs)
        new(value(lhs), value(rhs))
    end
end
Base.:(==)(a::Equation, b::Equation) = all(isequal.((a.lhs, a.rhs), (b.lhs, b.rhs)))
Base.hash(a::Equation, salt::UInt) = hash(a.lhs, hash(a.rhs, salt))

function Base.show(io::IO, eq::Equation)
    if hide_lhs(eq.lhs)
        show(io, eq.rhs)
    else
        print(io, eq.lhs, " ~ ", eq.rhs)
    end
end

SymbolicUtils.scalarize(eq::Equation) = scalarize(eq.lhs) .~ scalarize(eq.rhs)
SymbolicUtils.simplify(x::Equation; kw...) = simplify(x.lhs; kw...) ~ simplify(x.rhs; kw...)
# ambiguity
for T in [:Pair, :Any]
    @eval function SymbolicUtils.substitute(x::Equation, rules::$T; kw...)
        sub = substituter(rules)
        sub(x.lhs; kw...) ~ sub(x.rhs; kw...)
    end

    @eval function SymbolicUtils.substitute(eqs::Array{Equation}, rules::$T; kw...)
        sub = substituter(rules)
        sub.(lhss(eqs); kw...) .~ sub.(rhss(eqs); kw...)
    end
end

SymbolicUtils.substitute(nums::Array{Num}, rules; kw...) = substituter(rules).(nums; kw...)

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

julia> A ~ B
(broadcast(~, A, B))[1:3,1:3]

julia> A .~ 3x
(broadcast(~, A, 3x))[1:3,1:3]
```
"""
function Base.:~(lhs, rhs)
    if (isarraysymbolic(lhs) || isarraysymbolic(rhs)) && ((sl = size(lhs)) != (sr = size(rhs)))
        throw(ArgumentError("Cannot equate an array of different sizes. Got $sl and $sr."))
    else
        Equation(lhs, rhs)
    end
end
for T in [:Num, :Complex, :Number], S in [:Num, :Complex, :Number]
    (T != :Complex && S != :Complex) && continue
    @eval Base.:~(a::$T, b::$S) = let ar = value(real(a)), br = value(real(b)),
                                      ai = value(imag(a)), bi = value(imag(b))
        if ar isa Number && br isa Number && ai isa Number && bi isa Number
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

canonical_form(eq::Equation) = eq.lhs - eq.rhs ~ 0

function SymbolicUtils.search_variables!(buffer, eq::Equation; kw...)
    SymbolicUtils.search_variables!(buffer, eq.lhs; kw...)
    SymbolicUtils.search_variables!(buffer, eq.rhs; kw...)
end

function expand_derivatives(eq::Equation, simplify=false)
    return Equation(expand_derivatives(eq.lhs, simplify), expand_derivatives(eq.rhs, simplify))
end
