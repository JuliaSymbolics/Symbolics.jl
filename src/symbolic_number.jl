import SpecialFunctions: polygamma

# Atomic wrapper for symbolic scalar expressions whose symbolic type is numeric but not
# necessarily real. Unlike `Complex{Num}`, this wrapper does not decompose an expression
# into real and imaginary components. The wrapped `BasicSymbolic` remains one expression
# tree.
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

SymbolicNumber(x::SymbolicNumber) = x

SymbolicUtils.unwrap(x::SymbolicNumber) = x.val
SU.infer_vartype(::Type{SymbolicNumber}) = VartypeT
SymbolicUtils.symtype(x::SymbolicNumber) = symtype(unwrap(x))

# Route symbolic arithmetic through the raw BasicSymbolic algebra and only choose the
# wrapper after `promote_symtype` has determined the mathematical result domain. `//` is
# deliberately excluded: it constructs an exact Rational and is not generic division.
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
# `angle` is real-valued even when its argument is a general numeric symbolic scalar.
# SymbolicUtils does not currently register it among the standard monadic operations, so
# construct the real-typed symbolic call explicitly at this wrapper boundary.
Base.angle(x::SymbolicNumber) = Num(Term{VartypeT}(
    angle, ArgsT{VartypeT}((unwrap(x),)); type = Real, shape = SymbolicUtils.ShapeVecT()
))
Base.sincospi(x::SymbolicNumber) = (sinpi(x), cospi(x))
Base.numerator(x::SymbolicNumber) = wrap(numerator(unwrap(x)))
Base.denominator(x::SymbolicNumber) = wrap(denominator(unwrap(x)))
# `@number_methods` defines `^(::SymbolicNumber, ::Real)`, which intersects Base's
# integer/rational power methods. Keep those powers on the symbolic algebra explicitly.
Base.:^(x::SymbolicNumber, p::Integer) = wrap(unwrap(x)^p)
Base.:^(x::SymbolicNumber, p::Rational) = wrap(unwrap(x)^p)
# Base has a dedicated `ℯ ^ ::Number` method which intersects the generic symbolic
# exponent methods emitted above. Preserve the canonical exponential representation.
Base.:^(::Irrational{:ℯ}, x::SymbolicNumber) = wrap(exp(unwrap(x)))

# `polygamma(::Integer, ::Number)` in SpecialFunctions intersects the generic symbolic
# binary-function methods. This exact intersection keeps integer orders on the symbolic
# expression path without broadening the dispatch surface.
polygamma(m::Integer, x::SymbolicNumber) = wrap(polygamma(m, unwrap(x)))

# Base implements `cis(::Real)` through `sincos` followed by explicit `Complex`
# construction. That is appropriate for numerical values but would reintroduce Cartesian
# storage for `Num`. Canonically lower symbolic `cis` to the equivalent atomic scalar
# expression instead.
Base.cis(x::Num) = wrap(exp(im * unwrap(x)))

Base.iszero(x::SymbolicNumber) = SymbolicUtils._iszero(unwrap(x))
Base.isone(x::SymbolicNumber) = SymbolicUtils._isone(unwrap(x))
Base.zero(::SymbolicNumber) = SymbolicNumber(0)
Base.zero(::Type{SymbolicNumber}) = SymbolicNumber(0)
Base.one(::SymbolicNumber) = SymbolicNumber(1)
Base.one(::Type{SymbolicNumber}) = SymbolicNumber(1)

# `SymbolicNumber` is the wide numeric wrapper. Ordinary real values mixed with `Num`
# continue to promote to `Num` via `num.jl`; only the explicit Num/complex edge in
# `complex.jl` widens a real symbolic value to `SymbolicNumber`.
Base.promote_rule(::Type{SymbolicNumber}, ::Type{SymbolicNumber}) = SymbolicNumber
Base.promote_rule(::Type{T}, ::Type{SymbolicNumber}) where {T <: Number} = SymbolicNumber
Base.promote_rule(::Type{SymbolicNumber}, ::Type{T}) where {T <: Number} = SymbolicNumber
# Exact intersections with Base promotion rules keep Aqua ambiguity-free.
Base.promote_rule(::Type{Bool}, ::Type{SymbolicNumber}) = SymbolicNumber
Base.promote_rule(::Type{T}, ::Type{SymbolicNumber}) where {T <: AbstractIrrational} =
    SymbolicNumber
Base.promote_rule(::Type{Num}, ::Type{SymbolicNumber}) = SymbolicNumber
Base.promote_rule(::Type{SymbolicNumber}, ::Type{Num}) = SymbolicNumber
Base.convert(::Type{SymbolicNumber}, x::Number) = SymbolicNumber(x)

# Wrappers are representation boundaries, not distinct symbolic identities. Matching the
# wrapped expression's hash and `isequal` semantics is required by generic substitution,
# which recursively visits raw `BasicSymbolic` nodes while users naturally provide wrapped
# variables as dictionary keys.
Base.hash(x::SymbolicNumber, h::UInt) = hash(unwrap(x), h)::UInt
Base.isequal(a::SymbolicNumber, b::SymbolicNumber) = isequal(unwrap(a), unwrap(b))
Base.isequal(a::SymbolicNumber, b::BasicSymbolic) = isequal(unwrap(a), b)
Base.isequal(a::BasicSymbolic, b::SymbolicNumber) = isequal(a, unwrap(b))
Base.isequal(a::SymbolicNumber, b::Num) = isequal(unwrap(a), unwrap(b))
Base.isequal(a::Num, b::SymbolicNumber) = isequal(unwrap(a), unwrap(b))

function Base.show(io::IO, x::SymbolicNumber)
    warn_load_latexify()
    show(io, unwrap_const(unwrap(x)))
end

# Generic symbolic utilities must select the wrapper from the transformed expression's
# resulting symtype. In particular, a simplification is allowed to narrow a
# complex-capable expression to a provably real `Num`.
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
    SymbolicUtils.search_variables!(buffer, unwrap(expr); kw...)
end

SymbolicIndexingInterface.symbolic_type(::Type{SymbolicNumber}) = ScalarSymbolic()
SymbolicIndexingInterface.hasname(x::SymbolicNumber) = hasname(unwrap(x))
SymbolicIndexingInterface.getname(x::SymbolicNumber) = getname(unwrap(x))
function SymbolicIndexingInterface.symbolic_evaluate(x::SymbolicNumber, d::Dict; kw...)
    SymbolicIndexingInterface.symbolic_evaluate(unwrap(x), d; kw...)
end

function (s::SymbolicUtils.Substituter)(x::SymbolicNumber)
    wrap(s(unwrap(x)))
end
# The default substituter may widen a real `Num` expression when a replacement is complex.
# Specialize that concrete path so it selects the wrapper from the substituted symtype
# instead of forcing the result back through `Num`.
function (s::SymbolicUtils.DefaultSubstituter)(x::Num)
    wrap(s(unwrap(x)))
end

# High-level APIs historically reconstructed all symbolic derivatives as `Num`. Keep the
# existing real-valued paths unchanged and bridge only the wider scalar wrapper through the
# raw symbolic algorithms; `wrap` then recovers the mathematical result domain.
function derivative(O::SymbolicNumber, var; simplify = false, kwargs...)
    wrap(expand_derivatives(Differential(var)(unwrap(O)), simplify; kwargs...))
end
function derivative(O::AbstractArray{<:SymbolicNumber}, var; simplify = false, kwargs...)
    map(O) do o
        wrap(expand_derivatives(Differential(var)(unwrap(o)), simplify; kwargs...))
    end
end
derivative(f::Function, var::SymbolicNumber) = derivative(f(var), var)

function gradient(O::SymbolicNumber, vars::AbstractVector; simplify = false, kwargs...)
    map(vars) do var
        wrap(expand_derivatives(Differential(var)(unwrap(O)), simplify; kwargs...))
    end
end

function jacobian(
        ops::AbstractVector{<:SymbolicNumber}, vars::AbstractVector{<:SymbolicNumber};
        simplify = false, scalarize::Union{Val{true}, Val{false}} = Val(true), kwargs...
    )
    if scalarize isa Val{true}
        ops = Symbolics.scalarize(ops)
        vars = Symbolics.scalarize(vars)
    end
    raw_ops = unwrap.(ops)::Vector{SymbolicT}
    raw_vars = unwrap.(vars)::Vector{SymbolicT}
    return wrap.(jacobian(raw_ops, raw_vars; simplify, scalarize = Val(false), kwargs...))
end

function sparsejacobian_vals(
        ops::AbstractVector{<:SymbolicNumber}, vars::AbstractVector{<:SymbolicNumber},
        I::AbstractVector, J::AbstractVector; simplify::Bool = false, kwargs...
    )
    raw_ops = unwrap.(ops)::Vector{SymbolicT}
    raw_vars = unwrap.(vars)::Vector{SymbolicT}
    return [
        wrap(expand_derivatives(
            Differential(raw_vars[j])(raw_ops[i]), simplify; kwargs...
        )) for (i, j) in zip(I, J)
    ]
end

# `hessian(O, vars::Arr)` already materializes symbolic arrays with `collect(vars)`.
# Own the resulting concrete vector here; this avoids an ambiguity without depending on
# the later-defined `Arr` wrapper during package initialization.
function hessian(
        O::SymbolicNumber, vars::Vector{<:SymbolicNumber};
        simplify = false, kwargs...
    )
    return jacobian(
        gradient(O, vars; simplify, kwargs...), vars; simplify, kwargs...
    )
end

function sparsehessian(
        O::SymbolicNumber, vars::AbstractVector{<:SymbolicNumber};
        simplify::Bool = false, full::Bool = true, kwargs...
    )
    H = SparseArrays.sparse(hessian(O, vars; simplify, kwargs...))
    return full ? H : tril(H)
end

# Mirror the existing `Num`/`Complex{Num}` public LU bridge exactly. Constraining both the
# symbolic element type and the concrete/adjoint/transpose storage avoids intersections
# with LinearAlgebra's strided-matrix factorization methods.
function LinearAlgebra.lu(
        A::Union{
            Adjoint{<:SymbolicNumber}, Transpose{<:SymbolicNumber},
            Array{<:SymbolicNumber},
        }; check = true, kw...
    )
    sym_lu(A; check = check)
end

# Julia's numerical matrix exponential does not know about the atomic wrapper. Build the
# same symbolic matrix operation as the existing `Num`/`Complex{Num}` bridges and let
# `wrap` select the array wrapper from the resulting symbolic type at call time.
Base.exp(A::Matrix{SymbolicNumber}) = wrap(exp(SConst(A)))

# Linear expansion is wrapper-neutral internally; only the public scalar boundary needs to
# unwrap the variable and re-wrap the coefficient/remainder according to their symtypes.
function linear_expansion(t, x::SymbolicNumber)
    a, b, islinear = linear_expansion(t, unwrap(x))
    return wrap(a), wrap(b), islinear
end

# Match the existing `Num` convenience path for an equation array with a single scalar
# unknown. Without this bridge, dispatch falls into the scalar solver with an array-valued
# residual and eventually reaches `__solve(::Num, ::Arr{Equation,1}, ...)`.
function symbolic_linear_solve(
        eqs::AbstractArray, var::SymbolicNumber; simplify = false, check = true
    )
    return first(symbolic_linear_solve(eqs, [var]; simplify, check))
end

# A system whose unknowns are general numeric scalars may contain genuinely complex
# coefficients. Normalize equations to raw symbolic expressions before entering the
# existing array linear-expansion algorithm, which deliberately operates on SymbolicT.
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
