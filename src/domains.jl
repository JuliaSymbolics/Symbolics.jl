import DomainSets: Domain, Interval, AbstractInterval
import Symbolics: value, Sym, Term, Num

struct VarDomainPairing
    variables::Any
    domain::Domain
end

const DomainedVar = Union{Symbolic{<:Number}, Num}

Base.:∈(variable::DomainedVar, domain::Domain) = VarDomainPairing(value(variable), domain)
Base.:∈(variable::DomainedVar, domain::Interval) = VarDomainPairing(value(variable), domain)

# Construct Interval domain from a Tuple
function Base.:∈(variable::DomainedVar, domain::NTuple{2, Real})
    VarDomainPairing(variable, Interval(domain...))
end

# Multiple variables
function Base.:∈(variables::NTuple{N, DomainedVar}, domain::Domain) where {N}
    VarDomainPairing(value.(variables), domain)
end

function infimum(d::AbstractInterval{T}) where {T <: Num}
    leftendpoint(d)
end

function supremum(d::AbstractInterval{T}) where {T <: Num}
    rightendpoint(d)
end
