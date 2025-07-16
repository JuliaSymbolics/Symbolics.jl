@symbolic_wrap struct Num <: Real
    val
end

const RCNum = Union{Num, Complex{Num}}

unwrap(x::Num) = x.val

"""
    Num(val)

Wrap anything in a type that is a subtype of Real
"""
Num

const show_numwrap = Ref(false)

Num(x::Num) = x # ideally this should never be called
(n::Num)(args...) = Num(value(n)(map(value,args)...))
value(x) = unwrap(x)

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

Base.ifelse(x::Num,y,z) = Num(ifelse(value(x), value(y), value(z)))

Base.promote_rule(::Type{Bool}, ::Type{<:Num}) = Num
for C in [Complex, Complex{Bool}]
    @eval begin
        Base.:*(x::Num, z::$C) = Complex(x * real(z), x * imag(z))
        Base.:*(z::$C, x::Num) = Complex(real(z) * x, imag(z) * x)
        Base.:/(x::Num, z::$C) = let (a, b) = reim(z), den = a^2 + b^2
            Complex(x * a / den, -x * b / den)
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
Base.:^(::Irrational{:ℯ}, x::Num) = exp(x)

function Base.show(io::IO, z::Complex{<:Num})
    r, i = reim(z)
    compact = get(io, :compact, false)
    show(io, r)
    print(io, (compact ? "+" : " + ") * "(")
    show(io, i)
    print(io, ")*im")
end

# TODO: move this to SymbolicUtils
substitute(expr, s::Pair; kw...) = substituter([s[1] => s[2]])(expr; kw...)
substitute(expr, s::Vector; kw...) = substituter(s)(expr; kw...)

function _unwrap_callwithmeta(x)
    x = value(x)
    return x isa CallWithMetadata ? x.f : x
end
function subrules_to_dict(pairs)
    if pairs isa Pair
        pairs = (pairs,)
    end
    return Dict(_unwrap_callwithmeta(k) => value(v)  for (k, v) in pairs)
end
function substituter(pairs)
    dict = subrules_to_dict(pairs)
    (expr; kw...) -> SymbolicUtils.substitute(value(expr), dict; kw...)
end

SymbolicUtils.symtype(n::Num) = symtype(value(n))
Base.nameof(n::Num) = nameof(value(n))

Base.iszero(x::Num) = SymbolicUtils.fraction_iszero(unwrap(x))
Base.isone(x::Num) = SymbolicUtils.fraction_isone(unwrap(x))

import SymbolicUtils: <ₑ, Symbolic, Term, operation, arguments

Base.show(io::IO, n::Num) = show_numwrap[] ? print(io, :(Num($(value(n))))) : Base.show(io, value(n))

Base.promote_rule(::Type{<:Number}, ::Type{<:Num}) = Num
Base.promote_rule(::Type{BigFloat}, ::Type{<:Num}) = Num
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

Num(q::AbstractIrrational) = Num(Term(identity, [q]))

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

# Boolean operations
for (f, Domain) in [:(==) => :((AbstractFloat, Number)), :(!=) => :((AbstractFloat, Number)),
                    :(<=) => :((Real,)),   :(>=) => :((Real,)),
                    :(isless) => :((Real, AbstractFloat)),
                    :(<) => :((Real,)),   :(> ) => :((Real,)),
                    :(& )=> :((Bool,)),  :(| ) => :((Bool,)),
                    :xor => :((Bool,))]
    @eval @num_method Base.$f (val = $f(value(a), value(b)); val isa Bool ? val : Num(val)) $Domain
end

for f in [:!, :~]
    @eval Base.$f(x::Num) = (val = $f(value(x)); val isa Bool ? val : Num(val))
end
@num_method Base.isequal begin
  va = value(a)
  vb = value(b)
  if va isa SymbolicUtils.BasicSymbolic{Real} && vb isa SymbolicUtils.BasicSymbolic{Real}
    isequal(va, vb)::Bool
  else
    isequal(va, vb)::Bool
  end
end (AbstractFloat, Number, Symbolic)

Base.to_index(x::Num) = Base.to_index(value(x))

Base.hash(x::Num, h::UInt) = hash(value(x), h)::UInt

Base.convert(::Type{Num}, x::Symbolic{<:Number}) = Num(x)
Base.convert(::Type{Num}, x::Number) = Num(x)
Base.convert(::Type{Num}, x::Num) = x

Base.convert(::Type{T}, x::AbstractArray{Num}) where T <: Array{Num} = T(map(Num, x))
Base.convert(::Type{Sym}, x::Num) = value(x) isa Sym ? value(x) : error("cannot convert $x to Sym")

LinearAlgebra.lu(x::Union{Adjoint{<:RCNum},Transpose{<:RCNum},Array{<:RCNum}}; check=true, kw...) = sym_lu(x; check=check)

_iszero(x::Number) = iszero(x)
_isone(x::Number) = isone(x)
_iszero(::Symbolic) = false
_isone(::Symbolic) = false
_iszero(x::Num) = _iszero(value(x))::Bool
_isone(x::Num) = _isone(value(x))::Bool

Code.cse(x::Num) = Code.cse(unwrap(x))

## Documentation
# This method makes the docstring show all entries in the metadata dict associated with an instance of Num
function Base.Docs.getdoc(x::Num)
    x = unwrap(x)
    strings =
        ["A variable of type Symbolics.Num (Num wraps anything in a type that is a subtype of Real)";
        "# Metadata"]
    for (key, val) in collect(pairs(x.metadata))
        push!(strings, string(string(key), ": ", string(val)))
    end
    Markdown.parse(join(strings, "\n\n  "))
end

# https://github.com/JuliaSymbolics/Symbolics.jl/issues/1206#issuecomment-2271847091
"""
$(TYPEDSIGNATURES)

Return the alignment of printing `x` of type `Num`.

The alignment is a tuple `(left, right)` showing how many characters are needed 
on either side of an alignment feature. This function returns the text width
of `x` and `0` to avoid matching special characters, such as `e`and `f`, with 
the alignment algorithm in Julia Base, which leads to extra white spaces on the 
left of the characters when displaying array of symbolic variables.
"""
function Base.alignment(io::IO, x::Num)
    s = sprint(show, x, context = Base.nocolor(io), sizehint = 0)
    textwidth(s), 0
end
