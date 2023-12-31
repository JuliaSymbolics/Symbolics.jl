module PreallocationToolsSymbolicsExt

using PreallocationTools
import PreallocationTools: _restructure, get_tmp
using Symbolics, ForwardDiff

function get_tmp(dc::DiffCache, u::Type{X}) where {T,N, X<: ForwardDiff.Dual{T, Num, N}}
    if length(dc.du) > length(dc.any_du)
        resize!(dc.any_du, length(dc.du))
    end
    _restructure(dc.du, dc.any_du)
end

function get_tmp(dc::DiffCache, u::X) where {T,N, X<: ForwardDiff.Dual{T, Num, N}}
    if length(dc.du) > length(dc.any_du)
        resize!(dc.any_du, length(dc.du))
    end
    _restructure(dc.du, dc.any_du)
end

function get_tmp(dc::DiffCache, u::AbstractArray{X}) where {T,N, X<: ForwardDiff.Dual{T, Num, N}}
    if length(dc.du) > length(dc.any_du)
        resize!(dc.any_du, length(dc.du))
    end
    _restructure(dc.du, dc.any_du)
end

function get_tmp(dc::FixedSizeDiffCache, u::Type{X}) where {T,N, X<: ForwardDiff.Dual{T, Num, N}}
    if length(dc.du) > length(dc.any_du)
        resize!(dc.any_du, length(dc.du))
    end
    _restructure(dc.du, dc.any_du)
end

function get_tmp(dc::FixedSizeDiffCache, u::X) where {T,N, X<: ForwardDiff.Dual{T, Num, N}}
    if length(dc.du) > length(dc.any_du)
        resize!(dc.any_du, length(dc.du))
    end
    _restructure(dc.du, dc.any_du)
end

function get_tmp(dc::FixedSizeDiffCache, u::AbstractArray{X}) where {T,N, X<: ForwardDiff.Dual{T, Num, N}}
    if length(dc.du) > length(dc.any_du)
        resize!(dc.any_du, length(dc.du))
    end
    _restructure(dc.du, dc.any_du)
end

end
