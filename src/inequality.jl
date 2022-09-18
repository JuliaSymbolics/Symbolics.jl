"""
$(TYPEDEF)

An inequality relationship between two expressions.

# Fields
$(FIELDS)
"""
struct Inequality
    """The expression on the left-hand side of the inequality."""
    lhs
    """The expression on the right-hand side of the inequality."""
    rhs
    """The relational operator of the inequality."""
    relational_op
    function Inequality(lhs, rhs, relational_op)
        new(Symbolics.value(lhs), Symbolics.value(rhs), relational_op)
    end
end

Base.:(==)(a::Inequality, b::Inequality) = all([isequal(a.lhs, b.lhs), isequal(a.rhs, b.rhs), isequal(a.relational_op, b.relational_op)])
Base.hash(a::Inequality, salt::UInt) = hash(a.lhs, hash(a.rhs, hash(a.relational_op, salt)))

@enum RelationalOperator leq geq # strict less than or strict greater than are not supported by any solver

scalarize(ineq::Inequality) = Inequality(scalarize(ineq.lhs), scalarize(ineq.rhs), ineq.relational_op)

function Base.show(io::IO, ineq::Inequality)
    print(io, ineq.lhs, ineq.relational_op == leq ? " ≲ " : " ≳ ", ineq.rhs)
end

"""
$(TYPEDSIGNATURES)

Create an [`Inequality`](@ref) out of two [`Num`](@ref) instances, or an
`Num` and a `Number`.

# Examples

```jldoctest
julia> using Symbolics

julia> @variables x y;

julia> x ≲ y
x ≲ y

julia> x - y ≲ 0
x - y ≲ 0
```
"""
function ≲(lhs, rhs) 
    if isarraysymbolic(lhs) || isarraysymbolic(rhs)
        if isarraysymbolic(lhs) && isarraysymbolic(rhs)
            lhs .≲ rhs
        else
            throw(ArgumentError("Cannot relate an array with a scalar. Please use broadcast `.≲`."))
        end
    else
        Inequality(lhs, rhs, leq)
    end
end

"""
$(TYPEDSIGNATURES)

Create an [`Inequality`](@ref) out of two [`Num`](@ref) instances, or an
`Num` and a `Number`.

# Examples

```jldoctest
julia> using Symbolics

julia> @variables x y;

julia> x ≳ y
x ≳ y

julia> x - y ≳ 0
x - y ≳ 0
```
"""
function ≳(lhs, rhs)
    if isarraysymbolic(lhs) || isarraysymbolic(rhs)
        if isarraysymbolic(lhs) && isarraysymbolic(rhs)
            lhs .≳ rhs
        else
            throw(ArgumentError("Cannot relate an array with a scalar. Please use broadcast `.≳`."))
        end
    else
    Inequality(lhs, rhs, geq)
    end
end

function canonical_form(cs::Inequality; form=leq)
    # do we need to flip the operator?
    if cs.relational_op == form
        Inequality(cs.lhs - cs.rhs, 0, cs.relational_op)
    else
        Inequality(-(cs.lhs - cs.rhs), 0, cs.relational_op == leq ? geq : leq)
    end
end

SymbolicUtils.simplify(cs::Inequality; kw...) = 
    Inequality(simplify(cs.lhs; kw...), simplify(cs.rhs; kw...), cs.relational_op)
