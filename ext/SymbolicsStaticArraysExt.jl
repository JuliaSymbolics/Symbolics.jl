module SymbolicsStaticArraysExt

using Symbolics: Num, unwrap
import StaticArrays
import SymbolicUtils

function Base.:(/)(x::StaticArrays.SHermitianCompact{N, T}, y::Num) where {N, T <: Real}
    return SymbolicUtils.term(/, x, unwrap(y))
end
function Base.:(\)(x::Num, y::StaticArrays.SHermitianCompact{N, T}) where {N, T <: Real}
    return SymbolicUtils.term(\, unwrap(x), y)
end

end
