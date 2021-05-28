import DomainSets: Domain, Interval
import Symbolics: value, Sym, Term, Num

struct VarDomainPairing
  variables
  domain::Domain
end

Base.:∈(variable::Union{Sym,Term,Num},domain::Domain) = VarDomainPairing(value(variable),domain)
Base.:∈(variable::Union{Sym,Term,Num},domain::Interval) = VarDomainPairing(value(variable),domain)

# Construct Interval domain from a Tuple
Base.:∈(variable::Union{Sym,Term,Num},domain::NTuple{2,Real}) = VarDomainPairing(variable,Interval(domain...))

# Multiple variables
Base.:∈(variables::NTuple{N,Union{Sym,Term,Num}},domain::Domain) where N = VarDomainPairing(value.(variables),domain)
