module SymbolicsPreallocationToolsExt

using PreallocationTools
import PreallocationTools: get_tmp
using Symbolics, ForwardDiff

function _symbolic_tmp(dc::C, u::Type{X}) where {
        C <: Union{DiffCache, FixedSizeDiffCache}, X <: ForwardDiff.Dual,
    }
    return invoke(get_tmp, Tuple{C, Type{Y} where {Y <: Number}}, dc, u)
end

function get_tmp(dc::DiffCache, u::Type{X}) where {T, N, X <: ForwardDiff.Dual{T, Num, N}}
    return _symbolic_tmp(dc, u)
end

function get_tmp(dc::DiffCache, u::X) where {T, N, X <: ForwardDiff.Dual{T, Num, N}}
    return get_tmp(dc, X)
end

function get_tmp(dc::DiffCache, u::AbstractArray{X}) where {T, N, X <: ForwardDiff.Dual{T, Num, N}}
    return get_tmp(dc, X)
end

function get_tmp(dc::FixedSizeDiffCache, u::Type{X}) where {T, N, X <: ForwardDiff.Dual{T, Num, N}}
    return _symbolic_tmp(dc, u)
end

function get_tmp(dc::FixedSizeDiffCache, u::X) where {T, N, X <: ForwardDiff.Dual{T, Num, N}}
    return get_tmp(dc, X)
end

function get_tmp(dc::FixedSizeDiffCache, u::AbstractArray{X}) where {T, N, X <: ForwardDiff.Dual{T, Num, N}}
    return get_tmp(dc, X)
end

end
