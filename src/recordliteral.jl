"""
Maximum number of fields for which `@symstruct` generates symbolic constructor methods.
The number of methods is `2^nfields - 1`, so this bounds the cost of registering a struct.
Wider structs use [`record_literal`](@ref) explicitly.
"""
const RECORD_LITERAL_MAX_FIELDS = 8

"""
Argument types which mark a constructor call as symbolic, causing it to build a
[`RecordLiteral`](@ref) rather than calling the struct's own constructor.
"""
const RecordLiteralArg = Union{Num, Arr, SymStruct, SymbolicT}

"""
    $TYPEDEF

Operation used to represent a symbolic struct literal of type `T`. A term using this
operation carries symtype `T` and takes one argument per field of `T`, in `fieldnames(T)`
order.

This is the record counterpart of [`SymbolicUtils.array_literal`](@ref). An array whose
elements are symbolic is still an array, so `[p, q]` needs no special representation. A
struct is nominal and its fields are usually concretely typed, so `T(p, q)` cannot be
built as a value of `T` at all. A literal term stands in for that value: it names the
struct and its field expressions without ever calling `T`'s constructor, so fields may be
symbolic, and the constructor only runs during code generation, once the fields have
values.

Because the constructor is not called, any validation it performs is *not* applied to a
literal. Construction with fully concrete arguments is unaffected and still runs the real
constructor, including its validation.

Use [`record_literal`](@ref) to build these terms rather than constructing them directly.
"""
struct RecordLiteral{T} end

"""
    $TYPEDSIGNATURES

Return the struct type `f` is a literal of.
"""
record_type(::RecordLiteral{T}) where {T} = T

"""
    $TYPEDSIGNATURES

Return `true` if `x` is a symbolic struct literal (see [`RecordLiteral`](@ref)).
"""
is_record_literal(x) = false
function is_record_literal(x::SymbolicT)
    @match x begin
        BSImpl.Term(; f) => f isa RecordLiteral
        _ => false
    end
end

# A literal of fully concrete fields is just the value, so the operation is the constructor.
(::RecordLiteral{T})(args...) where {T} = T(args...)

SymbolicUtils.promote_symtype(::RecordLiteral{T}, args::SymbolicUtils.TypeT...) where {T} = T
function SymbolicUtils.promote_shape(::RecordLiteral, @nospecialize(args::SU.ShapeT...))
    return SU.ShapeVecT()
end

"""
    $TYPEDSIGNATURES

Build the symbolic struct literal of type `T` with field values `args`, given in
`fieldnames(T)` order. Returns a term of symtype `T`.

Prefer calling `T` directly - `@symstruct` registers constructor methods which produce a
literal when any argument is symbolic. This function is the explicit spelling, and is the
only one available for structs with more than `$(RECORD_LITERAL_MAX_FIELDS)` fields, where
generating those methods would be prohibitive.
"""
function record_literal(::Type{T}, args) where {T}
    nf = fieldcount(T)
    if length(args) != nf
        throw(ArgumentError(LazyString(
            "Cannot build a symbolic literal of `", T, "`: it has ", nf,
            " fields but ", length(args), " arguments were given.")))
    end
    cargs = ArgsT{VartypeT}()
    sizehint!(cargs, nf)
    for a in args
        push!(cargs, BSImpl.Const{VartypeT}(unwrap(a)))
    end
    return BSImpl.Term{VartypeT}(
        RecordLiteral{T}(), cargs; type = T, shape = SU.ShapeVecT())
end

record_literal(::Type{T}, args...) where {T} = record_literal(T, args)

function SymbolicUtils.show_call(io::IO, @nospecialize(f::RecordLiteral), x::SymbolicT)
    print(io, nameof(record_type(f)), "(")
    @match x begin
        BSImpl.Term(; args) => join(io, args, ", ")
    end
    print(io, ")")
end

function SymbolicUtils.Code.function_to_expr(
        @nospecialize(f::RecordLiteral), x::SymbolicT, st
    )
    out = get(st.rewrites, x, nothing)
    out === nothing || return out
    expr = Expr(:call, record_type(f))
    @match x begin
        BSImpl.Term(; args) => for a in args
            push!(expr.args, SymbolicUtils.Code.toexpr(a, st))
        end
    end
    return expr
end

"""
    $TYPEDSIGNATURES

The `codegen_function!` counterpart of the `function_to_expr` method above. A literal
lowers to an ordinary call to the struct's constructor, with each field's generated value
as an argument.
"""
function SymbolicUtils.Code.codegen_function!(
        @nospecialize(f::RecordLiteral), cs::SymbolicUtils.Code.CodegenState{T},
        expr::BasicSymbolic{T}, expr_idx::Integer
    ) where {T}
    result = Expr(:call, record_type(f))
    for arg in arguments(expr)
        push!(result.args, cs(arg))
    end
    return SymbolicUtils.Code.codegen!(cs, expr_idx, result)
end
