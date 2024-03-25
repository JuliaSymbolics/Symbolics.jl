import Base: Colon
# Note: we use Real because we don't yet have a
# specific wrapper for Integer expression types
# Real wrapper (Num) is the only one we have.
# Adding every new wrapper adds thousands of new
# methods to Base functions, we try to avoid this.
struct SymUnitRange <: AbstractUnitRange{Real}
    a
    b
end

Base.first(s::SymUnitRange) = s.a
Base.last(s::SymUnitRange) = s.b

(::Colon)(x::Num, y::Num)  = SymUnitRange(x, y)
(::Colon)(x::Num, y::Real) = SymUnitRange(x, y)
(::Colon)(x::Real, y::Num) = SymUnitRange(x, y)
Base.size(s::SymUnitRange) = (length(s),)
Base.length(s::SymUnitRange) = last(s) - (first(s) - 1)

struct Size
    size::Tuple
    unification_state::Dict
end


##
#
# @variables A[_1, _2] b[1:10] c[1:5]
#
# A * b --> infer that _2 = 1:10
# A * b + A * c --> fail.
# A * b + c' * A --> infer _1 = 1:5, _2 = 1:10
#
