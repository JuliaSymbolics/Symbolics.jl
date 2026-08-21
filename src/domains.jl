import DomainSets: Domain, Interval, AbstractInterval
import Symbolics: value, Sym, Term, Num

struct VarDomainPairing
  variables::BasicSymbolic{VartypeT}
  domain::Domain
end

const DomainedVar = Union{SymbolicT, Num}

Base.:∈(variable::DomainedVar,domain::Domain) = VarDomainPairing(unwrap(variable),domain)
Base.:∈(variable::DomainedVar,domain::Interval) = VarDomainPairing(unwrap(variable),domain)

# Construct Interval domain from a Tuple
Base.:∈(variable::DomainedVar,domain::NTuple{2,Real}) = VarDomainPairing(variable,Interval(domain...))

# Multiple variables
Base.:∈(variables::NTuple{N,DomainedVar},domain::Domain) where N = VarDomainPairing(unwrap.(variables),domain)

"""
    Symbolics.infimum(d::AbstractInterval)

Deprecated. Use `IntervalSets.infimum` (also re-exported by DomainSets) instead.

Return the left endpoint of the interval domain `d`.

This function was added so that an interval with `Num` endpoints, such as the domain of a
PDE variable written `x ∈ Interval(a, b)`, could be queried for its lower bound. That
turned out to be unnecessary: `IntervalSets.infimum` is generic over the endpoint type and
already handles `Num`, so the only effect of a separate `Symbolics.infimum` is to shadow it
whenever both Symbolics and DomainSets are in scope.

`Symbolics.infimum` now forwards to `IntervalSets.infimum` and warns; it will be removed in
the next breaking release. Callers should switch to

```julia
import DomainSets
infimum(d)
```

See also: [`Symbolics.supremum`](@ref).
"""
function infimum(d::AbstractInterval)
    Base.depwarn(
        "`Symbolics.infimum` is deprecated, use `IntervalSets.infimum` (re-exported by " *
            "DomainSets) instead.", :infimum
    )
    return IntervalSets.infimum(d)
end

"""
    Symbolics.supremum(d::AbstractInterval)

Deprecated. Use `IntervalSets.supremum` (also re-exported by DomainSets) instead.

Return the right endpoint of the interval domain `d`. This is the counterpart of
[`Symbolics.infimum`](@ref) and is deprecated for the same reason: `IntervalSets.supremum`
is generic over the endpoint type and already handles `Num`, so a separate Symbolics
function only shadows it.

`Symbolics.supremum` now forwards to `IntervalSets.supremum` and warns; it will be removed
in the next breaking release.
"""
function supremum(d::AbstractInterval)
    Base.depwarn(
        "`Symbolics.supremum` is deprecated, use `IntervalSets.supremum` (re-exported by " *
            "DomainSets) instead.", :supremum
    )
    return IntervalSets.supremum(d)
end
