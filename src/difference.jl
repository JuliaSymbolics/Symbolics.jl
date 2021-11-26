"""
$(TYPEDEF)

Represents a difference operator.

# Fields
$(FIELDS)

# Examples

```jldoctest
julia> using Symbolics

julia> @variables t;

julia> Δ = Difference(t; dt=0.01)
(::Difference) (generic function with 2 methods)
```
"""
struct Difference <: Operator
    """Fixed Difference"""
    t
    dt
    update::Bool
    Difference(t; dt, update=false) = new(value(t), dt, update)
end
(D::Difference)(t) = Term{symtype(t)}(D, [t])
(D::Difference)(t::Num) = Num(D(value(t)))
SymbolicUtils.promote_symtype(::Difference, t) = t
"""
$(SIGNATURES)

Represents a discrete update (shift) operator with the semantics
```
DiscreteUpdate(t; dt=0.01)(y) ~ y(t+dt)
```

# Examples

```jldoctest
julia> using Symbolics

julia> @variables t;

julia> U = DiscreteUpdate(t; dt=0.01)
(::Difference) (generic function with 2 methods)
```
"""
DiscreteUpdate(t; dt) = Difference(t; dt=dt, update=true)

Base.show(io::IO, D::Difference) = print(io, "Difference(", D.t, "; dt=", D.dt, ", update=", D.update, ")")

Base.:(==)(D1::Difference, D2::Difference) = isequal(D1.t, D2.t) && isequal(D1.dt, D2.dt) && isequal(D1.update, D2.update)
Base.hash(D::Difference, u::UInt) = hash(D.dt, hash(D.t, xor(u, 0x055640d6d952f101)))

Base.:^(D::Difference, n::Integer) = _repeat_apply(D, n)

"""
    hasdiff(O)

Returns true if the expression or equation `O` contains [`Difference`](@ref) terms (this include [`DiscreteUpdate`](@ref)).
"""
hasdiff(O) = recursive_hasoperator(Difference, O)