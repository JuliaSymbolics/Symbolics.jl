"""
$(TYPEDEF)

Represents a difference operator.

# Fields
$(FIELDS)

# Examples

```jldoctest
julia> using Symbolics

julia> @variables t;

julia> D = Difference(t; dt=0.01)
Difference(t; dt=0.01)
```
"""
struct Difference <: Function
    """Fixed Difference"""
    t
    dt
    Difference(t; dt) = new(value(t), dt)
end
(D::Difference)(t) = Term{symtype(t)}(D, [t])
(D::Difference)(t::Num) = Num(D(value(t)))
SymbolicUtils.promote_symtype(::Difference, t) = t

Base.show(io::IO, D::Difference) = print(io, "Difference(", D.t, "; dt=", D.dt, ")")

Base.:(==)(D1::Difference, D2::Difference) = isequal(D1.t, D2.t) && isequal(D1.dt, D2.dt)
Base.isequal(D1::Difference, D2::Difference) = isequal(D1.t, D2.t) && isequal(D1.dt, D2.dt)
Base.hash(D::Difference, u::UInt) = hash(D.dt, hash(D.t, xor(u, 0x055640d6d952f101)))
