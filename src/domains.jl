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

function infimum(d::AbstractInterval{<:Num})
    leftendpoint(d)
end

function supremum(d::AbstractInterval{<:Num})
    rightendpoint(d)
end
