import DomainSets: Domain, Interval
import Symbolics: value

struct VarDomainPairing
  variables
  domain::Domain
end

Base.:∈(variable::SymbolicUtils.Sym,domain::Domain) = VarDomainPairing(value(variable),domain)
Base.:∈(variable::SymbolicUtils.Sym,domain::Interval) = VarDomainPairing(value(variable),domain)

# Construct Interval domain from a Tuple
Base.:∈(variable::SymbolicUtils.Sym,domain::NTuple{2,Real}) = VarDomainPairing(variable,Interval(domain...))

# Multiple variables
Base.:∈(variables::NTuple{N,SymbolicUtils.Sym},domain::Domain) where N = VarDomainPairing(value.(variables),domain)
