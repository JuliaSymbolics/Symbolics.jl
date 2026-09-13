# Implement a few of the LogExpFunctions methods when those rely on boolean workflows.

LogExpFunctions.log1mexp(x::SymbolicScalar) = log(1 - exp(x))
LogExpFunctions.log1pexp(x::SymbolicScalar) = log(1 + exp(x))
LogExpFunctions.logexpm1(x::SymbolicScalar) = log(exp(x) - 1)
LogExpFunctions.logmxp1(x::SymbolicScalar) = log(x) - x + 1
for (f, op) in ((:logaddexp, +), (:logsubexp, -))
    @eval begin
        LogExpFunctions.$(f)(x::SymbolicScalar, y::Real) = log($(op)(exp(x), exp(y)))
        LogExpFunctions.$(f)(x::Real, y::SymbolicScalar) = log($(op)(exp(x), exp(y)))
        LogExpFunctions.$(f)(x::SymbolicScalar, y::SymbolicScalar) = log($(op)(exp(x), exp(y)))
    end
end
function LogExpFunctions.logsumexp(x::Union{AbstractVector{<:SymbolicScalar}, Arr})
    log(sum(exp, x; init = 0.0))
end
