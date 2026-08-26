module SymbolicsPreallocationToolsExt

using PreallocationTools
import PreallocationTools: get_tmp
using Symbolics, ForwardDiff

function _symbolic_tmp(dc::C,
        u::Union{X, Type{X}, AbstractArray{X}}) where {
        C <: Union{DiffCache, FixedSizeDiffCache}, X <: ForwardDiff.Dual
}
    return invoke(get_tmp, Tuple{C, Type{Y} where {Y <: Number}}, dc, X)
end

function get_tmp(dc::DiffCache,
        u::Union{X, Type{X}, AbstractArray{X}}) where {
        T, N, X <: ForwardDiff.Dual{T, Num, N}
}
    return _symbolic_tmp(dc, u)
end

function get_tmp(dc::FixedSizeDiffCache,
        u::Union{X, Type{X}, AbstractArray{X}}) where {
        T, N, X <: ForwardDiff.Dual{T, Num, N}
}
    return _symbolic_tmp(dc, u)
end

end
