"""
$(TYPEDEF)

The Fourier transform with respect to `t`, evaluated at the frequency `ω`, as an unevaluated
operator. `FourierTransform(t => ω)(f)` builds the term ``\\mathcal{F}[f](ω)`` as a
`Complex{Num}`, in the same way that [`Differential`](@ref) builds a derivative term;
[`expand_transforms`](@ref) applies the transform rules to it, and [`fourier`](@ref) does
both in one call.

The `convention` keyword selects the definition:

  - `:angular` (the default): ``F(ω) = ∫ f(t)\\, e^{-iωt}\\, dt``.
  - `:ordinary`: ``F(ν) = ∫ f(t)\\, e^{-2πiνt}\\, dt``, with ``ν`` in cycles per unit of `t`.
  - `:unitary`: ``F(ω) = (2π)^{-1/2} ∫ f(t)\\, e^{-iωt}\\, dt``.

`t` must be a real symbolic variable. The evaluation point `ω` must be real, free of `t`
and free of operator terms such as derivatives, integrals and other transforms.

The term is opaque to `substitute`, as `Differential` is to plain `substitute`: substitution
replaces a whole transform term but does not enter it. A derivative of a transform term stays
unevaluated.

# Fields
$(FIELDS)

# Examples

```julia
julia> @variables t ω;

julia> @syms f(t);

julia> FourierTransform(t => ω)(f(t))
FourierTransform(t => ω)(f(t))
```
"""
struct FourierTransform <: Operator
    """The integration (time) variable."""
    var::SymbolicT
    """The frequency at which the transform is evaluated."""
    point::SymbolicT
    """`:angular`, `:ordinary` or `:unitary`."""
    convention::Symbol
end

"""
$(TYPEDEF)

The inverse Fourier transform with respect to `ω`, evaluated at `t`, as an unevaluated
operator: `InverseFourierTransform(ω => t)(F)`. For each `convention` of
[`FourierTransform`](@ref) it is the inverse of the forward transform:

  - `:angular`: ``f(t) = (2π)^{-1} ∫ F(ω)\\, e^{iωt}\\, dω``.
  - `:ordinary`: ``f(t) = ∫ F(ν)\\, e^{2πiνt}\\, dν``.
  - `:unitary`: ``f(t) = (2π)^{-1/2} ∫ F(ω)\\, e^{iωt}\\, dω``.

It follows the same rules and restrictions as `FourierTransform`, with `ω` the integration
variable and `t` the evaluation point. [`inverse_fourier`](@ref) applies and expands it.

# Fields
$(FIELDS)
"""
struct InverseFourierTransform <: Operator
    """The integration (frequency) variable."""
    var::SymbolicT
    """The time at which the inverse transform is evaluated."""
    point::SymbolicT
    """`:angular`, `:ordinary` or `:unitary`."""
    convention::Symbol
end

const FourierOperator = Union{FourierTransform, InverseFourierTransform}

_isreal_symbolic(x::SymbolicT) = symtype(x) <: Real

function _fourier_operator(::Type{T}, (var, point)::Pair, convention) where {T}
    convention in (:angular, :ordinary, :unitary) ||
        throw(ArgumentError("Unknown Fourier convention `$(repr(convention))`. Use `:angular`, `:ordinary` or `:unitary`."))
    var = unwrap(var)
    var isa SymbolicT && issym(var) && _isreal_symbolic(var) ||
        throw(ArgumentError("The integration variable of a Fourier transform must be a real symbolic variable, got `$var`."))
    point = unwrap(point)
    point isa SymbolicT || (point = Const{VartypeT}(point))
    _isreal_symbolic(point) && !_depends_on(point, var) ||
        throw(ArgumentError("The evaluation point `$point` of a Fourier transform must be real, and free of `$var` and of operator terms."))
    return T(var, point, convention)
end
FourierTransform(p::Pair; convention = :angular) = _fourier_operator(FourierTransform, p, convention)
InverseFourierTransform(p::Pair; convention = :angular) = _fourier_operator(InverseFourierTransform, p, convention)

_transform_term(op::FourierOperator, x::SymbolicT) =
    Complex{Num}(Term{VartypeT}(op, SArgsT((x,)); type = Complex{Real}, shape = SymbolicUtils.ShapeVecT()))
(op::FourierOperator)(x::Num) = _transform_term(op, unwrap(x))
# the transform is linear over the complex numbers, so the integrand is always real-typed
(op::FourierOperator)(x::Complex{Num}) = _iszero(unwrap(imag(x))) ? op(real(x)) : op(real(x)) + im * op(imag(x))
(op::FourierOperator)(x::SymbolicT) = op(_transform_input(x))
SymbolicUtils.isbinop(::FourierOperator) = false

for T in (:FourierTransform, :InverseFourierTransform)
    @eval function Base.show(io::IO, op::$T)
        print(io, $(string(T)), "(", op.var, " => ", op.point)
        op.convention === :angular || print(io, "; convention = ", repr(op.convention))
        return print(io, ")")
    end
    @eval Base.isequal(a::$T, b::$T) = isequal(a.var, b.var) && isequal(a.point, b.point) && a.convention === b.convention
    @eval Base.:(==)(a::$T, b::$T) = isequal(a, b)
    @eval Base.hash(op::$T, h::UInt) = hash(op.convention, hash(op.point, hash(op.var, hash($T, h))))
end

is_fourier_term(x) = x isa SymbolicT && iscall(x) && operation(x) isa FourierOperator

# Whether `ex` may depend on `x`. Anything applied by an operation other than a function (a
# derivative, an integral, a transform) counts as depending on `x`, since its dependence is
# not tracked here.
function _depends_on(ex::SymbolicT, x::SymbolicT)
    isequal(ex, x) && return true
    iscall(ex) || return false
    f = operation(ex)
    f isa Function || f isa SymbolicT || return true
    return any(a -> _depends_on(a, x), arguments(ex))
end

# `(c, d)` with `u = c*x + d` and `c`, `d` free of `x`, or `nothing`. Only sums, products and
# quotients by factors free of `x` are accepted, so that `u` is affine by construction.
function _affine_parts(u::SymbolicT, x::SymbolicT)
    isequal(u, x) && return Num(1), Num(0)
    _depends_on(u, x) || return Num(0), Num(u)
    if isadd(u)
        parts = map(a -> _affine_parts(a, x), arguments(u))
        any(isnothing, parts) && return nothing
        return sum(first, parts), sum(last, parts)
    elseif ismul(u)
        args = arguments(u)
        i = findall(a -> _depends_on(a, x), args)
        length(i) == 1 || return nothing
        cd = _affine_parts(args[only(i)], x)
        cd === nothing && return nothing
        coeff = prod(Num, args[j] for j in eachindex(args) if j != only(i))
        return coeff * cd[1], coeff * cd[2]
    elseif isdiv(u)
        num, den = arguments(u)
        _depends_on(den, x) && return nothing
        cd = _affine_parts(num, x)
        cd === nothing && return nothing
        return cd[1] / Num(den), cd[2] / Num(den)
    end
    return nothing
end

# The transform is ∝ ∫ g(x) exp(iβkx) dx, evaluated at k.
function _kernel_beta(op::FourierOperator)
    β = op.convention === :ordinary ? -2 * Num(π) : Num(-1)
    return op isa FourierTransform ? β : -β
end
_at(op::T, point) where {T <: FourierOperator} = T(op.var, unwrap(point), op.convention)

function _transform(op::FourierOperator, g::SymbolicT)::Complex{Num}
    x = op.var
    _iszero(g) && return zero(Complex{Num})
    _depends_on(g, x) || return Num(g) * _transform_term(op, Const{VartypeT}(1))
    if isadd(g)
        return sum(a -> _transform(op, a), arguments(g))
    elseif ismul(g)
        args = arguments(g)
        coeff = filter(a -> !_depends_on(a, x), args)
        isempty(coeff) && return _transform_term(op, g)
        rest = filter(a -> _depends_on(a, x), args)
        return prod(Num, coeff) * _transform(op, unwrap(prod(Num, rest)))
    elseif isdiv(g)
        num, den = arguments(g)
        _depends_on(den, x) && return _transform_term(op, g)
        return _transform(op, num) / Num(den)
    elseif iscall(g) && operation(g) isa SymbolicT
        return _transform_affine(op, g)
    end
    return _transform_term(op, g)
end

# h(c x + d) ↔ e^{-iβkd/c} H(k/c) / |c|, for a symbolic function `h` of one argument
# depending on `x`, with `c` a finite nonzero real number and `d` real.
function _transform_affine(op::FourierOperator, g::SymbolicT)
    x = op.var
    args = arguments(g)
    i = findall(a -> _depends_on(a, x), args)
    length(i) == 1 || return _transform_term(op, g)
    cd = _affine_parts(args[only(i)], x)
    cd === nothing && return _transform_term(op, g)
    c, d = unwrap(cd[1]), unwrap(cd[2])
    SymbolicUtils.isconst(c) && unwrap_const(c) isa Real && isfinite(unwrap_const(c)) &&
        !iszero(unwrap_const(c)) && _isreal_symbolic(d) || return _transform_term(op, g)
    (isone(unwrap_const(c)) && _iszero(d)) && return _transform_term(op, g)
    ic = (c = unwrap_const(c)) isa Integer ? 1 // c : inv(c)
    newargs = collect(SymbolicT, args)
    newargs[only(i)] = x
    h = maketerm(SymbolicT, operation(g), newargs, metadata(g))
    k = Num(op.point)
    θ = -_kernel_beta(op) * k * Num(d) * ic
    phase = _iszero(unwrap(θ)) ? 1 : Complex{Num}(cos(θ), sin(θ))
    return phase * abs(ic) * _transform_term(_at(op, k * ic), h)
end

function _expand_transforms(ex::SymbolicT)::SymbolicT
    iscall(ex) || return ex
    f = operation(ex)
    args = arguments(ex)
    f isa FourierOperator && return ex
    if (f === real || f === imag) && is_fourier_term(only(args))
        z = only(args)
        return unwrap(f(_transform(operation(z), _expand_transforms(only(arguments(z))))))
    end
    newargs = map(_expand_transforms, args)
    all(i -> newargs[i] === args[i], eachindex(args)) && return ex
    return maketerm(SymbolicT, f, newargs, metadata(ex))
end

"""
    expand_transforms(ex)

Apply the Fourier transform rules to every [`FourierTransform`](@ref) and
[`Symbolics.InverseFourierTransform`](@ref) term in `ex`, innermost first. This is to
transforms what [`expand_derivatives`](@ref) is to `Differential`. The rules are:

  - linearity: sums, and factors that are provably free of the integration variable;
  - shift, scaling and reversal of a symbolic function, ``f(ct + d)``, with `c` a finite
    nonzero real number (not a symbol) and `d` real and free of `t`:
    ``e^{iωd/c} F(ω/c) / |c|`` in the angular convention.

A term no rule applies to stays unevaluated, so the result is exact. `ex` may be a `Num`, a
`Complex{Num}` or an array of them.
"""
expand_transforms(ex::Num) = Num(_expand_transforms(unwrap(ex)))
expand_transforms(ex::Complex{Num}) = Complex{Num}(expand_transforms(real(ex)), expand_transforms(imag(ex)))
expand_transforms(ex::AbstractArray) = map(expand_transforms, ex)
expand_transforms(ex) = ex

_transform_input(ex::Union{Num, Complex{Num}}) = ex
_transform_input(ex::Real) = Num(ex)
_transform_input(ex::Complex) = Complex{Num}(ex)
_transform_input(ex::SymbolicT) = symtype(ex) <: Complex ? Complex{Num}(ex) : Num(ex)

"""
    fourier(ex, t => ω; convention = :angular)

The Fourier transform of `ex` with respect to `t`, as a function of `ω`:
`expand_transforms(FourierTransform(t => ω; convention)(ex))`. The result is a
`Complex{Num}` in which any part no rule applies to stays an unevaluated
[`FourierTransform`](@ref) term. An array is transformed elementwise.

# Examples

```julia
julia> @variables t ω;

julia> @syms f(t);

julia> fourier(3f(t), t => ω)
3real(FourierTransform(t => ω)(f(t))) + 3im*imag(FourierTransform(t => ω)(f(t)))

julia> fourier(f(2t), t => ω)
(1//2)*real(FourierTransform(t => (1//2)*ω)(f(t))) + (1//2)*im*imag(FourierTransform(t => (1//2)*ω)(f(t)))
```
"""
function fourier(ex, p::Pair; convention = :angular)
    ex isa AbstractArray && return map(x -> fourier(x, p; convention), ex)
    return expand_transforms(FourierTransform(p; convention)(_transform_input(ex)))
end

"""
    inverse_fourier(F, ω => t; convention = :angular)

The inverse Fourier transform of `F` with respect to `ω`, as a function of `t`, under the
same `convention` as [`fourier`](@ref):
`expand_transforms(Symbolics.InverseFourierTransform(ω => t; convention)(F))`.
"""
function inverse_fourier(ex, p::Pair; convention = :angular)
    ex isa AbstractArray && return map(x -> inverse_fourier(x, p; convention), ex)
    return expand_transforms(InverseFourierTransform(p; convention)(_transform_input(ex)))
end
