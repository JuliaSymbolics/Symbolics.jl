"""
$(TYPEDEF)

An inequality relationship between two expressions.

# Fields
$(FIELDS)
"""
struct Inequality
    """The expression on the left-hand side of the inequality."""
    lhs::BasicSymbolic{VartypeT}
    """The expression on the right-hand side of the inequality."""
    rhs::BasicSymbolic{VartypeT}
    """The relational operator of the inequality, [`leq`](@ref) or [`geq`](@ref)."""
    relational_op
    function Inequality(lhs, rhs, relational_op)
        new(unwrap(lhs), unwrap(rhs), relational_op)
    end
end

Base.:(==)(a::Inequality, b::Inequality) = isequal(a.lhs, b.lhs) && isequal(a.rhs, b.rhs) && isequal(a.relational_op, b.relational_op)::Bool
Base.hash(a::Inequality, salt::UInt) = hash(a.lhs, hash(a.rhs, hash(a.relational_op, salt)))

@enum RelationalOperator leq geq # strict less than or strict greater than are not supported by any solver

"""
    leq

The relational operator of `lhs ≲ rhs`: `lhs` is less than or equal to `rhs`. Pass it
as the third argument of [`Inequality`](@ref), or compare it with the `relational_op`
field of one.
"""
leq

"""
    geq

The relational operator of `lhs ≳ rhs`: `lhs` is greater than or equal to `rhs`. Pass
it as the third argument of [`Inequality`](@ref), or compare it with the `relational_op`
field of one.
"""
geq

is_array_operand(x) = SU.is_array_shape(SU.shape(x))

function check_inequality_shapes(lhs, rhs, op)
    (is_array_operand(lhs) && is_array_operand(rhs)) || return nothing
    (SU.shape(lhs) isa Unknown || SU.shape(rhs) isa Unknown) && return nothing
    if (sl = size(lhs)) != (sr = size(rhs))
        sym = op == leq ? "≲" : "≳"
        throw(ArgumentError("Cannot relate arrays of different sizes $sl and $sr with \
            `$sym`. Use broadcast `.$sym` for elementwise inequalities."))
    end
    return nothing
end

function SymbolicUtils.scalarize(ineq::Inequality)
    lhs = scalarize(ineq.lhs)
    rhs = scalarize(ineq.rhs)
    if lhs isa AbstractArray || rhs isa AbstractArray
        check_inequality_shapes(lhs, rhs, ineq.relational_op)
        Inequality.(lhs, rhs, Ref(ineq.relational_op))
    else
        Inequality(lhs, rhs, ineq.relational_op)
    end
end

function Base.show(io::IO, ineq::Inequality)
    warn_load_latexify()
    print(io, ineq.lhs, ineq.relational_op == leq ? " ≲ " : " ≳ ", ineq.rhs)
end

"""
$(TYPEDSIGNATURES)

Create an [`Inequality`](@ref) out of two [`Num`](@ref) instances, or an `Num` and a `Number`.
Unicode `≲` can be typed by writing `\\lesssim` then pressing tab in the Julia REPL, and in many editors.

Either side may be an array. Two arrays must have the same size, and a scalar side bounds
every element of the array side. The result is a single array-valued `Inequality`; use
broadcast `.≲` for an array of scalar inequalities instead.

# Examples

```jldoctest
julia> using Symbolics

julia> @variables x y z[1:2];

julia> x ≲ y
x ≲ y

julia> x - y ≲ 0
x - y ≲ 0

julia> z ≲ x
z ≲ x

julia> scalarize(z ≲ x)
2-element Vector{Inequality}:
 z[1] ≲ x
 z[2] ≲ x
```
"""
function ≲(lhs, rhs)
    check_inequality_shapes(lhs, rhs, leq)
    Inequality(lhs, rhs, leq)
end

"""
$(TYPEDSIGNATURES)

Create an [`Inequality`](@ref) out of two [`Num`](@ref) instances, or an `Num` and a `Number`.
Unicode `≳` can be typed by writing `\\gtrsim` then pressing tab in the Julia REPL, and in many editors.

Either side may be an array, as for [`≲`](@ref).

# Examples

```jldoctest
julia> using Symbolics

julia> @variables x y z[1:2] w[1:2];

julia> x ≳ y
x ≳ y

julia> x - y ≳ 0
x - y ≳ 0

julia> z ≳ w
z ≳ w
```
"""
function ≳(lhs, rhs)
    check_inequality_shapes(lhs, rhs, geq)
    Inequality(lhs, rhs, geq)
end

inequality_difference(a, b) = is_array_operand(a) || is_array_operand(b) ? a .- b : a - b

function canonical_form(cs::Inequality; form=leq)
    # do we need to flip the operator?
    if cs.relational_op == form
        Inequality(inequality_difference(cs.lhs, cs.rhs), 0, cs.relational_op)
    else
        Inequality(inequality_difference(cs.rhs, cs.lhs), 0, cs.relational_op == leq ? geq : leq)
    end
end

function SymbolicUtils.search_variables!(buffer, ineq::Inequality; kw...)
    search_variables!(buffer, ineq.lhs; kw...)
    search_variables!(buffer, ineq.rhs; kw...)
end

SymbolicUtils.simplify(cs::Inequality; kw...) = 
    Inequality(simplify(cs.lhs; kw...), simplify(cs.rhs; kw...), cs.relational_op)

function (s::SymbolicUtils.Substituter)(x::Inequality)
    Inequality(s(x.lhs), s(x.rhs), x.relational_op)
end
