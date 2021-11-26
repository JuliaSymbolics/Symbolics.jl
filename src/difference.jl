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

See also [`Shift`](@ref)
"""
struct Difference <: Function
    """Fixed Difference"""
    t
    dt
    shift::Bool
    Difference(t; dt, shift=false) = new(value(t), dt, shift)
end
(D::Difference)(t) = Term{symtype(t)}(D, [t])
(D::Difference)(t::Num) = Num(D(value(t)))
SymbolicUtils.promote_symtype(::Difference, t) = t
"""
$(SIGNATURES)

Represents a discrete time-shift operator with the semantics
```
Shift(t; dt=0.01)(y) ~ y(t+dt)
```

# Examples

```jldoctest
julia> using Symbolics

julia> @variables t;

julia> U = Shift(t; dt=0.01)
(::Difference) (generic function with 2 methods)
```
"""
Shift(t; dt) = Difference(t; dt=dt, shift=true)

Base.@deprecate DiscreteUpdate(t; dt) Shift(t; dt)

Base.show(io::IO, D::Difference) = print(io, "Difference(", D.t, "; dt=", D.dt, ", shift=", D.shift, ")")

Base.:(==)(D1::Difference, D2::Difference) = isequal(D1.t, D2.t) && isequal(D1.dt, D2.dt) && isequal(D1.shift, D2.shift)
Base.hash(D::Difference, u::UInt) = hash(D.dt, hash(D.t, xor(u, 0x055640d6d952f101)))

Base.:^(D::Difference, n::Integer) = _repeat_apply(D, n)