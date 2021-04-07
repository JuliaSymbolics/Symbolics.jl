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

SciMLBase.issymbollike(::Num) = true
SciMLBase.issymbollike(::SymbolicUtils.Symbolic) = true

SymbolicUtils.@number_methods(
                              Num,
                              Num(f(value(a))),
                              Num(f(value(a), value(b))),
                              [conj, real, transpose]
                             )
Base.conj(x::Num) = x
Base.transpose(x::Num) = x

Base.eps(::Type{Num}) = Num(0)
Base.typemin(::Type{Num}) = Num(-Inf)
Base.typemax(::Type{Num}) = Num(Inf)
Base.float(x::Num) = x

IfElse.ifelse(x::Num,y,z) = Num(IfElse.ifelse(value(x), value(y), value(z)))

Base.promote_rule(::Type{Bool}, ::Type{<:Num}) = Num
for C in [Complex, Complex{Bool}]
    @eval begin
        Base.:*(x::Num, z::$C) = Complex(x * real(z), x * imag(z))
        Base.:*(z::$C, x::Num) = Complex(real(z) * x, imag(z) * x)
        Base.:/(x::Num, z::$C) = let (a, b) = reim(z), den = a^2 + b^2
            Complex(x * a / den, x * b / den)
        end
        Base.:/(z::$C, x::Num) = Complex(real(z) / x, imag(z) / x)
        Base.:+(x::Num, z::$C) = Complex(x + real(z), imag(z))
        Base.:+(z::$C, x::Num) = Complex(real(z) + x, imag(z))
        Base.:-(x::Num, z::$C) = Complex(x - real(z), -imag(z))
        Base.:-(z::$C, x::Num) = Complex(real(z) - x, imag(z))
    end
end

function Base.inv(z::Complex{Num})
    a, b = reim(z)
    den = a^2 + b^2
    Complex(a/den, -b/den)
end
function Base.:/(x::Complex{Num}, y::Complex{Num})
    a, b = reim(x)
    c, d = reim(y)
    den = c^2 + d^2
    Complex((a*c + b*d)/den, (b*c - a*d)/den)
end
Base.:^(z::Complex{Num}, n::Integer) = Base.power_by_squaring(z, n)

function Base.show(io::IO, z::Complex{<:Num})
    r, i = reim(z)
    compact = get(io, :compact, false)
    show(io, r)
    print(io, (compact ? "+" : " + ") * "(")
    show(io, i)
    print(io, ")*im")
end

SymbolicUtils.simplify(n::Num; kw...) = Num(SymbolicUtils.simplify(value(n); kw...))
SymbolicUtils.expand(n::Num) = Num(SymbolicUtils.expand(value(n)))
"""
    substitute(expr, s::Dict)

Performs the substitution on `expr` according to rule(s) `s`.
# Examples
```julia
julia> @parameters t
(t,)
julia> @variables x y z(t)
(x, y, z(t))
julia> ex = x + y + sin(z)
(x + y) + sin(z(t))
julia> substitute(ex, Dict([x => z, sin(z) => z^2]))
(z(t) + y) + (z(t) ^ 2)
```
"""
substitute(expr::Num, s::Pair; kw...) = Num(substituter(s)(value(expr); kw...)) # backward compat
substitute(expr::Num, s::Vector; kw...) = Num(substituter(s)(value(expr); kw...))
substitute(expr::Num, s::Dict; kw...) = Num(substituter(s)(value(expr); kw...))
# TODO: move this to SymbolicUtils
substitute(expr, s::Pair; kw...) = substituter([s[1] => s[2]])(expr; kw...)
substitute(expr, s::Vector; kw...) = substituter(s)(expr; kw...)

substituter(pair::Pair) = substituter((pair,))
function substituter(pairs)
    dict = Dict(value(k) => value(v)  for (k, v) in pairs)
    (expr; kw...) -> SymbolicUtils.substitute(expr, dict; kw...)
end

SymbolicUtils.symtype(n::Num) = symtype(n.val)

function Base.iszero(x::Num)
    x = value(x)
    x isa Number && iszero(x) && return true
    _x = SymbolicUtils.to_mpoly(x)[1]
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

@num_method Base.isless  (val = isless(tosymbolic(a), tosymbolic(b)); val isa Bool ? val : Num(val)) (Real,)
@num_method Base.:(<)    (val = tosymbolic(a) < tosymbolic(b)       ; val isa Bool ? val : Num(val)) (Real,)
@num_method Base.:(<=)   (val = tosymbolic(a) <= tosymbolic(b)      ; val isa Bool ? val : Num(val)) (Real,)
@num_method Base.:(>)    (val = tosymbolic(a) > tosymbolic(b)       ; val isa Bool ? val : Num(val)) (Real,)
@num_method Base.:(>=)   (val = tosymbolic(a) >= tosymbolic(b)      ; val isa Bool ? val : Num(val)) (Real,)
@num_method Base.:(==)   (val = tosymbolic(a) == tosymbolic(b)      ; val isa Bool ? val : Num(val)) (AbstractFloat,Number)
@num_method Base.isequal isequal(tosymbolic(a), tosymbolic(b)) (AbstractFloat, Number, Symbolic)

Base.hash(x::Num, h::UInt) = hash(value(x), h)

Base.convert(::Type{Num}, x::Symbolic{<:Number}) = Num(x)
Base.convert(::Type{Num}, x::Number) = Num(x)
Base.convert(::Type{Num}, x::Num) = x

Base.convert(::Type{<:Array{Num}}, x::AbstractArray) = map(Num, x)
Base.convert(::Type{<:Array{Num}}, x::AbstractArray{Num}) = x
Base.convert(::Type{Sym}, x::Num) = value(x) isa Sym ? value(x) : error("cannot convert $x to Sym")

LinearAlgebra.lu(x::Union{Adjoint{<:Num},Transpose{<:Num},Array{<:Num}}; check=true, kw...) = sym_lu(x; check=check)

_iszero(x::Number) = iszero(x)
_isone(x::Number) = isone(x)
_iszero(::Symbolic) = false
_isone(::Symbolic) = false
_iszero(x::Num) = _iszero(value(x))
_isone(x::Num) = _isone(value(x))

SymbolicUtils.Code.toexpr(x::Num) = SymbolicUtils.Code.toexpr(value(x))

SymbolicUtils.setmetadata(x::Num, t, v) = Num(SymbolicUtils.setmetadata(value(x), t, v))
SymbolicUtils.getmetadata(x::Num, t) = SymbolicUtils.getmetadata(value(x), t)
SymbolicUtils.hasmetadata(x::Num, t) = SymbolicUtils.hasmetadata(value(x), t)
