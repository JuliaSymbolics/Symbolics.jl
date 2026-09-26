"""
    NAMESPACE_SEPARATOR

Character used to separate nested symbolic names when Symbolics displays a
namespace-qualified variable.

The separator is `₊` and is used by variable construction and display code.
"""
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
    lhs::BasicSymbolic{VartypeT}
    """The expression on the right-hand side of the equation."""
    rhs::BasicSymbolic{VartypeT}
    function Equation(lhs, rhs)
        new(unwrap(lhs), unwrap(rhs))
    end
end
Base.:(==)(a::Equation, b::Equation) = isequal(a.lhs, b.lhs) && isequal(a.rhs, b.rhs)
Base.hash(a::Equation, salt::UInt) = hash(a.lhs, hash(a.rhs, salt))

function Base.show(io::IO, eq::Equation)
    warn_load_latexify()
    if hide_lhs(unwrap_const(eq.lhs))::Bool
        show(io, unwrap_const(eq.rhs))
    else
        print(io, eq.lhs, " ~ ", eq.rhs)
    end
end

SymbolicUtils.scalarize(eq::Equation, args...) = scalarize(eq.lhs, args...) .~ scalarize(eq.rhs, args...)
SymbolicUtils.simplify(x::Equation; kw...) = simplify(x.lhs; kw...) ~ simplify(x.rhs; kw...)
function (s::SymbolicUtils.Substituter)(eq::Equation)
    s(eq.lhs) ~ s(eq.rhs)
end

lhss(xs) = map(x->x.lhs, xs)
rhss(xs) = map(x->x.rhs, xs)

"""
$(TYPEDEF)

The real-part and imaginary-part [`Equation`](@ref)s of a complex-valued
equation, as produced by [`split_complex_equation`](@ref).

Unlike the plain `Vector{Equation}` returned by a complex [`~`](@ref), a
`SplitComplexEquation` records that its two equations came from splitting one
complex equation, so it can be told apart from a user-written pair of real
equations with `x isa SplitComplexEquation`. `eqs.original` recovers the
unsplit complex equation.

It is a read-only two-element `AbstractVector{Equation}` with `eqs[1] ==
eqs.real_eq` and `eqs[2] == eqs.imag_eq`: iteration, indexing, `length`,
`vcat` and `reduce(vcat, ...)` work, but `push!` and `setindex!` do not. Call
`collect(eqs)` where a `Vector{Equation}` is required.

# Fields
$(FIELDS)
"""
struct SplitComplexEquation <: AbstractVector{Equation}
    """The equation between the real parts of the two sides."""
    real_eq::Equation
    """The equation between the imaginary parts of the two sides."""
    imag_eq::Equation
    """The unsplit equation between the original complex-valued sides."""
    original::Equation
end

Base.size(::SplitComplexEquation) = (2,)

function Base.getindex(eqs::SplitComplexEquation, i::Int)
    i == 1 && return eqs.real_eq
    i == 2 && return eqs.imag_eq
    throw(BoundsError(eqs, i))
end

Base.iterate(eqs::SplitComplexEquation, state::Int = 1) =
    state > 2 ? nothing : (eqs[state], state + 1)

"""
    split_complex_equation(lhs, rhs)

Equate `lhs` and `rhs` like `lhs ~ rhs`, but return a
[`SplitComplexEquation`](@ref) instead of a `Vector{Equation}` when the
equation splits into real-part and imaginary-part equations.

When only the real parts or only the imaginary parts contain symbols, the
result is a single [`Equation`](@ref), exactly as for `~`. Throws an error when
neither side contains symbols.

# Examples

```jldoctest
julia> using Symbolics

julia> @variables x y;

julia> eqs = split_complex_equation(x + im * y, 1 + 2im)
2-element SplitComplexEquation:
 x ~ 1
 y ~ 2

julia> eqs isa SplitComplexEquation
true

julia> split_complex_equation(x, 1 + 0im)
x ~ 1
```
"""
function split_complex_equation(a::Number, b::Number)
    ar, br = value(real(a)), value(real(b))
    ai, bi = value(imag(a)), value(imag(b))
    return if ar isa Number && br isa Number && ai isa Number && bi isa Number
        error("Equation $a ~ $b does not contain any symbols")
    elseif ar isa Number && br isa Number
        ai ~ bi
    elseif ai isa Number && bi isa Number
        ar ~ br
    else
        SplitComplexEquation(ar ~ br, ai ~ bi, Equation(a, b))
    end
end

"""
$(TYPEDSIGNATURES)

Create an [`Equation`](@ref) out of two [`Num`](@ref) instances, or an
`Num` and a `Number`.

When either side is complex-valued and both the real and imaginary parts
contain symbols, `~` returns a `Vector{Equation}` holding the real-part and
imaginary-part equations. Use [`split_complex_equation`](@ref) to get the same
pair marked as a [`SplitComplexEquation`](@ref).

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
A ~ B

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
    @eval function Base.:~(a::$T, b::$S)
        eq = split_complex_equation(a, b)
        return eq isa SplitComplexEquation ? Equation[eq.real_eq, eq.imag_eq] : eq
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
