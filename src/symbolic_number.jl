import SpecialFunctions: polygamma

@symbolic_wrap struct SymbolicNumber <: Number
    val::BasicSymbolic{VartypeT}

    function SymbolicNumber(ex::BasicSymbolic{VartypeT})
        @assert symtype(ex) <: Number
        return new(Const{VartypeT}(ex))
    end

    function SymbolicNumber(ex::Number)
        return new(Const{VartypeT}(unwrap(ex)))
    end
end

"""
    SymbolicNumber(ex)

Wrap a symbolic scalar whose numeric domain is not known to be real, such as a variable
declared with `@variables z::Complex` or the expression `im * x`. A `SymbolicNumber` is a
`Number`, while [`Num`](@ref) is a `Real`.

The wrapped expression is atomic: it is not split into real and imaginary parts.
Operations build the raw symbolic expression and choose the wrapper from its symtype, so a
result known to be real, such as `real(z)`, `imag(z)` or `abs(z)`, is a `Num`, and any
other numeric result is a `SymbolicNumber`.

`Complex{Num}` remains the explicit Cartesian representation. Mixing it with a
`SymbolicNumber` promotes to `SymbolicNumber`, and `Complex(reim(z)...)` converts a
`SymbolicNumber` `z` to Cartesian form. See the manual page on complex numbers.
"""
SymbolicNumber

SymbolicNumber(x::SymbolicNumber) = x

SymbolicUtils.unwrap(x::SymbolicNumber) = x.val
SU.infer_vartype(::Type{SymbolicNumber}) = VartypeT
SymbolicUtils.symtype(x::SymbolicNumber) = symtype(unwrap(x))

# `//` constructs an exact `Rational`; it is not symbolic division.
SymbolicUtils.@number_methods(
    SymbolicNumber,
    wrap(f(unwrap(a))),
    wrap(f(unwrap(a), unwrap(b))),
    [conj, real, imag, transpose, //],
)

Base.conj(x::SymbolicNumber) = wrap(conj(unwrap(x)))
Base.real(x::SymbolicNumber) = wrap(real(unwrap(x)))
Base.imag(x::SymbolicNumber) = wrap(imag(unwrap(x)))
Base.transpose(x::SymbolicNumber) = wrap(transpose(unwrap(x)))
Base.adjoint(x::SymbolicNumber) = wrap(adjoint(unwrap(x)))
# SymbolicUtils has no real-valued `angle` rule, so the term is built with `type = Real`.
Base.angle(x::SymbolicNumber) = Num(
    Term{VartypeT}(
        angle, ArgsT{VartypeT}((unwrap(x),)); type = Real, shape = SymbolicUtils.ShapeVecT()
    )
)
Base.sincospi(x::SymbolicNumber) = (sinpi(x), cospi(x))
Base.numerator(x::SymbolicNumber) = wrap(numerator(unwrap(x)))
Base.denominator(x::SymbolicNumber) = wrap(denominator(unwrap(x)))
# Resolve ambiguities with Base's integer and rational powers.
Base.:^(x::SymbolicNumber, p::Integer) = wrap(unwrap(x)^p)
Base.:^(x::SymbolicNumber, p::Rational) = wrap(unwrap(x)^p)
# Resolve an ambiguity with Base's `ℯ^x`.
Base.:^(::Irrational{:ℯ}, x::SymbolicNumber) = wrap(exp(unwrap(x)))

# Resolve an ambiguity with `polygamma(::Integer, ::Number)` from SpecialFunctions.
polygamma(m::Integer, x::SymbolicNumber) = wrap(polygamma(m, unwrap(x)))

# Base builds `cis` from `sincos`, which would split the result into Cartesian parts.
Base.cis(x::Num) = wrap(exp(im * unwrap(x)))

Base.iszero(x::SymbolicNumber) = SymbolicUtils._iszero(unwrap(x))
Base.isone(x::SymbolicNumber) = SymbolicUtils._isone(unwrap(x))
Base.zero(::SymbolicNumber) = SymbolicNumber(0)
Base.zero(::Type{SymbolicNumber}) = SymbolicNumber(0)
Base.one(::SymbolicNumber) = SymbolicNumber(1)
Base.one(::Type{SymbolicNumber}) = SymbolicNumber(1)

Base.promote_rule(::Type{SymbolicNumber}, ::Type{SymbolicNumber}) = SymbolicNumber
Base.promote_rule(::Type{T}, ::Type{SymbolicNumber}) where {T <: Number} = SymbolicNumber
Base.promote_rule(::Type{SymbolicNumber}, ::Type{T}) where {T <: Number} = SymbolicNumber
# Resolve ambiguities with Base promotion rules.
Base.promote_rule(::Type{Bool}, ::Type{SymbolicNumber}) = SymbolicNumber
Base.promote_rule(::Type{T}, ::Type{SymbolicNumber}) where {T <: AbstractIrrational} =
    SymbolicNumber
Base.promote_rule(::Type{Num}, ::Type{SymbolicNumber}) = SymbolicNumber
Base.promote_rule(::Type{SymbolicNumber}, ::Type{Num}) = SymbolicNumber
Base.convert(::Type{SymbolicNumber}, x::Number) = SymbolicNumber(x)

# Substitution looks up wrapped keys while visiting raw nodes.
Base.hash(x::SymbolicNumber, h::UInt) = hash(unwrap(x), h)::UInt
Base.isequal(a::SymbolicNumber, b::SymbolicNumber) = isequal(unwrap(a), unwrap(b))
Base.isequal(a::SymbolicNumber, b::BasicSymbolic) = isequal(unwrap(a), b)
Base.isequal(a::BasicSymbolic, b::SymbolicNumber) = isequal(a, unwrap(b))
Base.isequal(a::SymbolicNumber, b::Num) = isequal(unwrap(a), unwrap(b))
Base.isequal(a::Num, b::SymbolicNumber) = isequal(unwrap(a), unwrap(b))

function Base.show(io::IO, x::SymbolicNumber)
    warn_load_latexify()
    return show(io, unwrap_const(unwrap(x)))
end

SymbolicUtils.simplify(x::SymbolicNumber; kw...) = wrap(SymbolicUtils.simplify(unwrap(x); kw...))
SymbolicUtils.simplify_fractions(x::SymbolicNumber; kw...) = wrap(SymbolicUtils.simplify_fractions(unwrap(x); kw...))
SymbolicUtils.expand(x::SymbolicNumber) = wrap(SymbolicUtils.expand(unwrap(x)))
SymbolicUtils.Code.toexpr(x::SymbolicNumber) = SymbolicUtils.Code.toexpr(unwrap(x))
SymbolicUtils.setmetadata(x::SymbolicNumber, t, v) = wrap(SymbolicUtils.setmetadata(unwrap(x), t, v))
SymbolicUtils.getmetadata(x::SymbolicNumber, t) = SymbolicUtils.getmetadata(unwrap(x), t)
SymbolicUtils.hasmetadata(x::SymbolicNumber, t) = SymbolicUtils.hasmetadata(unwrap(x), t)
Broadcast.broadcastable(x::SymbolicNumber) = x
SymbolicUtils.scalarize(x::SymbolicNumber) = wrap(SymbolicUtils.scalarize(unwrap(x)))

function SymbolicUtils.search_variables!(buffer, expr::SymbolicNumber; kw...)
    return SymbolicUtils.search_variables!(buffer, unwrap(expr); kw...)
end

SymbolicIndexingInterface.symbolic_type(::Type{SymbolicNumber}) = ScalarSymbolic()
SymbolicIndexingInterface.hasname(x::SymbolicNumber) = hasname(unwrap(x))
SymbolicIndexingInterface.getname(x::SymbolicNumber) = getname(unwrap(x))
function SymbolicIndexingInterface.symbolic_evaluate(x::SymbolicNumber, d::Dict; kw...)
    return SymbolicIndexingInterface.symbolic_evaluate(unwrap(x), d; kw...)
end

function (s::SymbolicUtils.Substituter)(x::SymbolicNumber)
    return wrap(s(unwrap(x)))
end

# Restricting the storage types avoids ambiguities with LinearAlgebra's strided `lu`.
function LinearAlgebra.lu(
        A::Union{
            Adjoint{<:SymbolicNumber}, Transpose{<:SymbolicNumber},
            Array{<:SymbolicNumber},
        }; check = true, kw...
    )
    return sym_lu(A; check = check)
end

Base.exp(A::Matrix{SymbolicNumber}) = wrap(exp(SConst(A)))

function linear_expansion(t, x::SymbolicNumber)
    a, b, islinear = linear_expansion(t, unwrap(x))
    return wrap(a), wrap(b), islinear
end

function symbolic_linear_solve(
        eqs::AbstractArray, var::SymbolicNumber; simplify = false, check = true
    )
    return first(symbolic_linear_solve(eqs, [var]; simplify, check))
end

function symbolic_linear_solve(
        eqs::AbstractArray, vars::AbstractArray{<:SymbolicNumber};
        simplify = false, check = true
    )
    raw_eqs = SymbolicT[
        eq isa Equation ? unwrap(eq.rhs) - unwrap(eq.lhs) : unwrap(eq)
            for eq in eqs
    ]
    A, b, islinear = linear_expansion(raw_eqs, unwrap.(vars))
    check && @assert islinear
    islinear || return nothing

    Aw = SymbolicNumber.(A)
    rhs = SymbolicNumber.(-b)
    sol = copy(rhs)
    LinearAlgebra.ldiv!(sym_lu(Aw), sol)
    return simplify ? SymbolicUtils.simplify_fractions.(sol) : sol
end
