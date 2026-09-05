"""
    $TYPEDEF

Operation used to represent a lazy symbolic zero of symtype `T`. Terms using this operation
take no arguments and carry symtype `T`, so that differentiating an expression which does
not depend on the differentiation variable produces a zero of the *same* symtype as the
expression, rather than a scalar `0`.

This exists because not every type has a materializable zero value. A `SymStruct` type may
have array fields of unknown shape, and calling the struct's constructor to build a concrete
zero instance is not safe in general. A lazy zero sidesteps both problems: it never commits
to a size, and never invokes a user-defined constructor. Field access and indexing fold
through it via [`symbolic_zero`](@ref).

Use [`symbolic_zero`](@ref) to construct zeros rather than building these terms directly;
it returns a constant where one exists, and only falls back to this representation when
necessary.
"""
struct SymbolicZero{T} end

"""
    $TYPEDSIGNATURES

Return the symtype that `f` represents the zero of.
"""
symzero_type(::SymbolicZero{T}) where {T} = T

SymbolicUtils.promote_symtype(::SymbolicZero{T}) where {T} = T

function (::SymbolicZero{T})() where {T}
    return symbolic_zero_term(T, shape_from_type(T))
end

"""
    $TYPEDSIGNATURES

Return `true` if `x` is a lazy symbolic zero term (see [`SymbolicZero`](@ref)). Note that
this is `false` for zeros which are represented as constants, so it is not a general
"is this zero" predicate.
"""
is_symbolic_zero(x) = false
function is_symbolic_zero(x::SymbolicT)
    @match x begin
        BSImpl.Term(; f) => f isa SymbolicZero
        _ => false
    end
end

"""
    $TYPEDSIGNATURES

Construct the lazy zero term of symtype `T` with shape `sh`.
"""
function symbolic_zero_term(::Type{T}, @nospecialize(sh)) where {T}
    return BSImpl.Term{VartypeT}(
        SymbolicZero{T}(), ArgsT{VartypeT}(); type = T, shape = sh)
end

"""
    $TYPEDSIGNATURES

Return `true` if a symbolic zero of symtype `T` can be constructed by
[`symbolic_zero`](@ref).

Numbers are always zeroable. An array type is zeroable if its element type is, and a
symbolic struct type is zeroable if every one of its fields is. Everything else - notably
`String` fields, which have no additive identity - is not.
"""
is_zeroable(::Type{T}) where {T} = _is_zeroable(T, Base.IdSet{Any}())

function _is_zeroable(@nospecialize(T::Type), seen::Base.IdSet{Any})
    T <: Number && return true
    # Guards against types which are recursive through an array field, e.g.
    # `struct Tree; kids::Vector{Tree}; end`.
    T in seen && return true
    push!(seen, T)
    T <: AbstractArray && return _is_zeroable(eltype(T), seen)
    is_symstruct_type(T) || return false
    for fT in fieldtypes(T)
        _is_zeroable(fT, seen) || return false
    end
    return true
end

"""
    $TYPEDSIGNATURES

Return the symbolic zero of symtype `T` with shape `sh`.

Where a zero can be represented as a symbolic constant it is, so that it folds and
simplifies like any other constant:

- Number types give `zero(T)`.
- Array types with a fully known shape and numeric elements give a concrete array of zeros,
  which folds to a scalar zero when indexed.

Otherwise the zero is represented lazily as a [`SymbolicZero`](@ref) term of symtype `T`.
This covers symbolic structs, whose constructors must not be called speculatively, and
array types whose shape is not known. Accessing a field of a lazy zero struct recursively
yields the zero of that field, so the laziness is not observable through field access.

Throws an `ArgumentError` if `T` is not zeroable; use [`is_zeroable`](@ref) to check first.
"""
function symbolic_zero(::Type{T}, @nospecialize(sh) = shape_from_type(T)) where {T}
    if T <: Number
        return BSImpl.Const{VartypeT}(zero(T))
    elseif T <: AbstractArray
        eT = eltype(T)
        if eT <: Number && sh isa SU.ShapeVecT && all(ax -> first(ax) == 1, sh)
            val = zeros(eT, map(length, sh)...)
            # Only usable if it reproduces `T` exactly - `T` may be a static array or some
            # other `AbstractArray` that `zeros` does not construct.
            typeof(val) === T && return BSImpl.Const{VartypeT}(val)
        end
        is_zeroable(eT) || _throw_not_zeroable(T)
        return symbolic_zero_term(T, sh)
    elseif is_symstruct_type(T)
        is_zeroable(T) || _throw_not_zeroable(T)
        return symbolic_zero_term(T, sh)
    end
    _throw_not_zeroable(T)
end

function _throw_not_zeroable(@nospecialize(T::Type))
    throw(ArgumentError("""
        Cannot construct a symbolic zero of type `$T`, since it has no additive identity. \
        This typically happens when differentiating a symbolic struct with a field whose \
        type is neither a `Number`, an array, nor another symbolic struct."""))
end

function SymbolicUtils.show_call(io::IO, @nospecialize(f::SymbolicZero), x::SymbolicT)
    print(io, "zero(", symzero_type(f), ")")
end

function SymbolicUtils.Code.function_to_expr(@nospecialize(f::SymbolicZero), x::SymbolicT, st)
    out = get(st.rewrites, x, nothing)
    out === nothing || return out
    throw(ArgumentError("""
        Cannot generate code for the symbolic zero of type `$(symzero_type(f))`, since it \
        has no concrete value. Access its fields or index into it to obtain values which \
        can be code generated."""))
end
