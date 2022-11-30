import DomainSets: Domain, Interval, AbstractInterval
import Symbolics: value, Sym, Term, Num

struct VarDomainPairing
  variables
  domain::Domain
end

for D in [:Domain, :Interval, :AbstractInterval]
    @eval Base.:∈(variable::Union{Sym,Term,Num},domain::$D) = VarDomainPairing(value(variable),domain)
end

# Construct Interval domain from a Tuple
Base.:∈(variable::Union{Sym,Term,Num},domain::NTuple{2,Real}) = VarDomainPairing(variable,Interval(domain...))

# Multiple variables
Base.:∈(variables::NTuple{N,Union{Sym,Term,Num}},domain::Domain) where N = VarDomainPairing(value.(variables),domain)

function infimum(d::AbstractInterval{T}) where T <: Num
    leftendpoint(d)
end

function supremum(d::AbstractInterval{T}) where T <: Num
    rightendpoint(d)
end
