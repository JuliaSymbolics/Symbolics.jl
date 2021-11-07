function optimize(x; kws...)
    # if iswrapped(x)
        return wrap(SymbolicUtils.optimize(unwrap(x); kws...))
    # end
    # return SymbolicUtils.optimize(x)
end

function optimize(x::AbstractArray; kws...)
    return wrap.(SymbolicUtils.optimize(unwrap.(x)))
end

import Base.map 
Base.map(f::typeof(optimize), x::AbstractArray) = optimize(x)

function getcost(x)
    SymbolicUtils.getcost(unwrap(x))
end