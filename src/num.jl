"""
    Num(val)

Wrap anything in a type that is a subtype of Real
"""
struct Num <: Real
    val
end

const show_numwrap = Ref(false)

Num(x::Num) = x # ideally this should never be called
(n::Num)(args...) = Num(value(n)(map(value,args)...))
value(x) = x
value(x::Num) = x.val


using SymbolicUtils: to_symbolic
SymbolicUtils.to_symbolic(n::Num) = value(n)
SymbolicUtils.@number_methods(Num,
                              Num(f(value(a))),
                              Num(f(value(a), value(b))))

SymbolicUtils.simplify(n::Num; kw...) = Num(SymbolicUtils.simplify(value(n); kw...))

SymbolicUtils.symtype(n::Num) = symtype(n.val)

function Base.iszero(x::Num)
    _x = SymbolicUtils.to_mpoly(value(x))[1]
    return (_x isa Number || _x isa SymbolicUtils.MPoly) && iszero(_x)
end

import SymbolicUtils: <ₑ, Symbolic, Term, operation, arguments

Base.show(io::IO, n::Num) = show_numwrap[] ? print(io, :(Num($(value(n))))) : Base.show(io, value(n))

Base.promote_rule(::Type{<:Number}, ::Type{<:Num}) = Num
Base.promote_rule(::Type{<:Symbolic{<:Number}}, ::Type{<:Num}) = Num
function Base.getproperty(t::Union{Add, Mul, Pow, Term}, f::Symbol)
    if f === :op
        Base.depwarn("`x.op` is deprecated, use `operation(x)` instead", :getproperty, force=true)
        operation(t)
    elseif f === :args
        Base.depwarn("`x.args` is deprecated, use `arguments(x)` instead", :getproperty, force=true)
        arguments(t)
    else
        getfield(t, f)
    end
end
<ₑ(s::Num, x) = value(s) <ₑ value(x)
<ₑ(s, x::Num) = value(s) <ₑ value(x)
<ₑ(s::Num, x::Num) = value(s) <ₑ value(x)

for T in (Integer, Rational)
    @eval Base.:(^)(n::Num, i::$T) = Num(value(n)^i)
end

macro num_method(f, expr, Ts=nothing)
    if Ts === nothing
        Ts = [Any]
    else
        @assert Ts.head == :tuple
        # e.g. a tuple or vector
        Ts = Ts.args
    end

    ms = [quote
              $f(a::$T, b::$Num) = $expr
              $f(a::$Num, b::$T) = $expr
          end for T in Ts]
    quote
        $f(a::$Num, b::$Num) = $expr
        $(ms...)
    end |> esc
end

"""
    tosymbolic(a::Union{Sym,Num}) -> Sym{Real}
    tosymbolic(a::T) -> T
"""
tosymbolic(a::Num) = tosymbolic(value(a))
tosymbolic(a::Sym) = Sym{symtype(a)}(nameof(a)) # unwrap stuff like Parameter{<:Number}
tosymbolic(a) = a

@num_method Base.isless isless(tosymbolic(a), tosymbolic(b)) (Real,)
@num_method Base.:(<) (tosymbolic(a) < tosymbolic(b)) (Real,)
@num_method Base.:(<=) (tosymbolic(a) <= tosymbolic(b)) (Real,)
@num_method Base.:(>) (tosymbolic(a) > tosymbolic(b)) (Real,)
@num_method Base.:(>=) (tosymbolic(a) >= tosymbolic(b)) (Real,)
@num_method Base.isequal isequal(tosymbolic(a), tosymbolic(b)) (AbstractFloat, Number, Symbolic)
@num_method Base.:(==) tosymbolic(a) == tosymbolic(b) (AbstractFloat,Number)

Base.hash(x::Num, h::UInt) = hash(value(x), h)

Base.convert(::Type{Num}, x::Symbolic{<:Number}) = Num(x)
Base.convert(::Type{Num}, x::Number) = Num(x)
Base.convert(::Type{Num}, x::Num) = x

Base.convert(::Type{<:Array{Num}}, x::AbstractArray) = map(Num, x)
Base.convert(::Type{<:Array{Num}}, x::AbstractArray{Num}) = x
Base.convert(::Type{Sym}, x::Num) = value(x) isa Sym ? value(x) : error("cannot convert $x to Sym")

LinearAlgebra.lu(x::Array{Num}; check=true, kw...) = sym_lu(x; check=check)

