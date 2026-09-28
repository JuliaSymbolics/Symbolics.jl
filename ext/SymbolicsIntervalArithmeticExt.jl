module SymbolicsIntervalArithmeticExt

using IntervalArithmetic
using Symbolics
using Symbolics: Num
import SymbolicUtils

const NumTypes = IntervalArithmetic.NumTypes

# Resolve the promote_rule ambiguity between
#   IntervalArithmetic: promote_rule(::Type{Interval{T}}, ::Type{S}) where {S<:Real}
#   Symbolics:          promote_rule(::Type{T}, ::Type{Num}) where {T<:Number}
# Prefer Num so intervals can appear as constant coefficients in symbolic expressions.
# See https://github.com/JuliaSymbolics/Symbolics.jl/issues/1157
Base.promote_rule(::Type{Interval{T}}, ::Type{Num}) where {T <: NumTypes} = Num
Base.promote_rule(::Type{Num}, ::Type{Interval{T}}) where {T <: NumTypes} = Num

# IntervalArithmetic purposely throws on inconclusive `==` / `isone` / `iszero`.
# SymbolicUtils uses those predicates (and Dict-keyed caches via `isequal`) when
# building products, so map them to the interval-safe equivalents.
Base.isone(x::Interval) = IntervalArithmetic.isthinone(x)
Base.iszero(x::Interval) = IntervalArithmetic.isthinzero(x)
Base.isequal(x::Interval, y::Interval) = IntervalArithmetic.isequal_interval(x, y)

# Printing of `a * x` checks `_isunit(a)` / `_isminus(a)` via `a == ±1`, which also throws.
SymbolicUtils._isunit(x::Interval) = IntervalArithmetic.isthinone(x)
SymbolicUtils._isminus(x::Interval) = IntervalArithmetic.isthin(x, -1)

end
