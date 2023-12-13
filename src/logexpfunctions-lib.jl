# Implement a few of the LogExpFunctions methods when those rely on boolean workflows.

LogExpFunctions.log1mexp(x::RCNum) = log(1 - exp(x))
LogExpFunctions.log1pexp(x::RCNum) = log(1 + exp(x))
LogExpFunctions.logexpm1(x::RCNum) = log(exp(x) - 1)
LogExpFunctions.logmxp1(x::RCNum) = log(x) - x + 1
for (f, op) in ((:logaddexp, +), (:logsubexp, -))
    @eval begin
        LogExpFunctions.$(f)(x::RCNum, y::Real) = log($(op)(exp(x), exp(y)))
        LogExpFunctions.$(f)(x::Real, y::RCNum) = log($(op)(exp(x), exp(y)))
        LogExpFunctions.$(f)(x::RCNum, y::RCNum) = log($(op)(exp(x), exp(y)))
    end
end
function LogExpFunctions.logsumexp(x::Union{AbstractVector{<:Num}, Arr})
    log(sum(exp, x; init = 0.0))
end
