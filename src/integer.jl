# using Symbolics: wrap, unwrap, operation, symtype, istree, arguments, Symbolic, getname, value, set_where, @symbolic_wrap
# using Symbolics

@symbolic_wrap struct IntegerTerm <: Integer
    n
end

IntegerTerm(x::IntegerTerm) = x # ideally this should never be called
(n::IntegerTerm)(args...) = IntegerTerm(value(n)(map(value,args)...))
value(x) = x
value(x::IntegerTerm) = unwrap(x)

Symbolics.unwrap(x::IntegerTerm) = x.n
Symbolics.symtype(a::IntegerTerm) = Integer
Symbolics.istree(a::IntegerTerm) = true
Symbolics.operation(a::IntegerTerm) = Integer
Symbolics.arguments(a::IntegerTerm) = [a.n]


Base.promote_rule(::Type{<:Integer}, ::Type{<:IntegerTerm}) = IntegerTerm
Base.promote_rule(::Type{<:Symbolic{<:Integer}}, ::Type{<:IntegerTerm}) = IntegerTerm
Base.promote_rule(::Type{<:Number}, ::Type{<:IntegerTerm}) = Num
Base.promote_rule(::Type{<:Num}, ::Type{<:IntegerTerm}) = Num


# SymbolicUtils.@number_methods(
#                               IntegerTerm,
#                               unwrap(f(wrap(a))),
#                               unwrap(f(wrap(a), wrap(b))),
#                              )

SymbolicUtils.@number_methods(
                              IntegerTerm,
                              IntegerTerm(f(value(a))),
                              IntegerTerm(f(value(a), value(b))),
                              [conj, real, transpose]
                             )