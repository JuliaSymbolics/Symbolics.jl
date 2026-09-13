"""
    $TYPEDEF

Wrapper type for symbolic structs. Requires that the wrapped struct type `T` be registered
with [`@symstruct`](@ref). After registration, `@variables` can be used to create the
symbolic struct.

```julia
# Here, `record` has type `SymStruct{Record}`
@variables record::Record
```

`getproperty` access on this is a symbolic operation, and returns an expression performing
the appropriate field access. This can only wrap concrete struct types (`isconcretetype(T)`
must be `true`). `getproperty` on this struct leverages `fieldnames` and `fieldtypes`.
Thus, will thus not respect custom `getproperty` methods on the wrapped struct type.

`SymStruct{T}` subtypes `AbstractVector{Symbolics.SymbolicT}`. For a fully registered type
where all fields (including those of nested symbolic structs) are of known size, `SymStruct{T}`
can be linearly indexed to access fields in order. Array fields are indexed linearly.
"""
struct SymStruct{T}
    sym::SymbolicT

    function SymStruct{T}(x::SymbolicT) where {T}
        @assert isconcretetype(T)
        # Ensure wrapped types have been registered as such
        @assert wrapper_type(T) === SymStruct{T}
        # Ensure that the symbolic represents the correct type
        @assert symtype(x) === T
        new{T}(x)
    end
end

is_wrapper_type(::Type{SymStruct}) = false
is_wrapper_type(::Type{S}) where {T, S <: SymStruct{<:T}} = true
wraps_type(::Type{S}) where {T, S <: SymStruct{T}} = T
iswrapped(::SymStruct{T}) where {T} = true

"""
    issymstruct(x) -> Bool

Return `true` if `x` is a `SymStruct` or a symbolic expression whose `symtype` has been
registered via [`@symstruct`](@ref). Return `false` otherwise.
"""
issymstruct(x) = false
issymstruct(x::SymStruct) = true
function issymstruct(x::SymbolicT)
    T = symtype(x)
    return has_symwrapper(T)::Bool && wrapper_type(T)::DataType === SymStruct{T}
end

SymbolicUtils.unwrap(x::SymStruct) = getfield(x, 1)
SymbolicUtils.infer_vartype(::Type{SymStruct{T}}) where {T} = VartypeT

function Base.isequal(a::SymStruct{T}, b::SymStruct{T}) where {T}
    isequal(unwrap(a), unwrap(b))
end
Base.isequal(a::SymStruct{T}, b::SymbolicT) where {T} = isequal(unwrap(a), b)
Base.isequal(a::SymbolicT, b::SymStruct{T}) where {T} = isequal(a, unwrap(b))
Base.hash(a::SymStruct{T}, h::UInt) where {T} = hash(unwrap(a), h)

function field_shape end

"""
    @symstruct Foo{T1, T2, ...}
    @symstruct Foo{T1, T2, ...} begin
      # options...
    end

A macro which enables using type `Foo` with `SymStruct` as a symbolic struct. The first
argument to the macro must be the struct type, with all type parameters named. The optional
second argument is an optional `begin..end` block containing options that influence the
behavior of the macro. The following options are allowed:

- `shape(:field) = # expression`. For array fields, the shape of the field cannot be
  inferred from the type. In case the type of the field can be inferred from the
  type, it can be specified using this syntax. The expression must evaluate to an object of type
  `Union{SymbolicUtils.Unknown, AbstractVector{UnitRange{Int}}, Tuple{Vararg{UnitRange{Int}}}}`.
  The expression has access to the concrete type of the struct being accessed, with all
  type parameters available as declared in the first argument.

For example, given the following struct:

```julia
struct Record{T}
  x::Int
  y::Real
  z::T
end
```

It can be registered as

```julia
# Note: the type parameter must be declared, but the name itself does not matter
@symstruct Record{V} begin
# If `V` is an `AbstractVector` then field `z` is a 3-vector. Otherwise, it is a scalar.
  shape(:z) = V <: AbstractVector ? [1:3] : ()
end
```

Now,

```julia
@variables rec::Record{Int} rec2::Record{Vector{Int}}
```

`rec.x`, `rec2.x` will be `Num`s with symtype `Int`. `rec.y` and `rec2.y` will be `Num`s
with symtype `Real`. `rec.z` will be a `Num` with symtype `Int`. `rec2.z` will be an
`Arr{Num, 1}` with symtype `Vector{Int}` and shape `[1:2]`.

In case the shape of a field is not provided, it will be inferred from the type. For
`AbstractArray` subtypes, it will be `SymbolicUtils.Unknown(ndims(arr_type))`. Otherwise,
it will be treated as a scalar.
"""
macro symstruct(T, opts = Expr(:block))
    block = Expr(:block)
    where_args = Expr[]
    nocurly_name = T
    if Meta.isexpr(T, :curly)
        for x in @view(T.args[2:end])
            push!(where_args, esc(x))
        end
        nocurly_name = T.args[1]
    end
    T = esc(T)
    nocurly_name = esc(nocurly_name)
    temp_typevar = :S
    push!(block.args, quote
        function (::$(typeof(has_symwrapper)))(::Type{$temp_typevar}) where {$(where_args...), $temp_typevar <: $T}
            true
        end
        function (::$(typeof(wrapper_type)))(::Type{$temp_typevar}) where {$(where_args...), $temp_typevar <: $T}
            isconcretetype($temp_typevar) ? $SymStruct{$temp_typevar} : $SymStruct{<:$temp_typevar}
        end
    end)

    @assert Meta.isexpr(opts, :block) """
    Options to `@symstruct` must be specified as a `begin...end` block. Got $opts.
    """
    for stmt in opts.args
        stmt isa LineNumberNode && continue
        @assert Meta.isexpr(stmt, :(=)) """
        Each option to `@symstruct` must be of the form `option(args...) = value`. \
        Got $stmt.
        """
        head, val = stmt.args
        @assert Meta.isexpr(head, :call) """
        Each option to `@symstruct` must be of the form `option(args...) = value`. \
        Got $head instead of `option(args...)`.
        """
        opt = head.args[1]
        args = @view(head.args[2:end])
        if opt === :shape
            @assert length(args) == 1 """
            The `shape` option must be of the form `shape(:field_name) = value`. Instead \
            of a single argument `:field_name`, multiple arguments $args were found.
            """
            @assert args[1] isa QuoteNode """
            The field name provided to the `shape` option must be a literal `Symbol`.
            Found `$(args[1])`.
            """
            field = args[1]
            push!(block.args, __field_shape_expr(T, field, where_args, val))
        else
            error("Unsupported option $opt.")
        end
    end

    return block
end

function __field_shape_expr(T::Union{Symbol, Expr}, field::QuoteNode,
                            where_args::Vector{Expr}, val::Union{Expr, Symbol})
    quote
        function (::$(typeof(field_shape)))(sym::Type{S}, ::Val{$field}) where {$(where_args...), S <: $T}
            val = $(esc(val))
            if val isa $(SymbolicUtils.Unknown)
                return val
            elseif val isa $(SymbolicUtils.ShapeVecT)
                return val
            elseif val isa $(AbstractVector{UnitRange{Int}})
                return $(SymbolicUtils.ShapeVecT)(val)
            elseif val isa $(Tuple{Vararg{UnitRange{Int}}})
                return $(SymbolicUtils.ShapeVecT)(val)
            else
                error("""
                Invalid usage of `@symstruct` macro for type $($T). The shape for field \
                $($field) was specified incorrectly. The result of the expression must be \
                one of `SymbolicUtils.Unknown`, `AbstractVector{UnitRange{Int}}` or \
                `Tuple{Vararg{UnitRange{Int}}}`. Found a value of type $(typeof(val)).
                """)
            end
        end
    end
end

# Generated `if..elseif..else` chain for `getproperty`.
@generated function Base.getproperty(sym::SymStruct{T}, name::Symbol) where {T}
    chain = Expr(:if)
    cur = chain
    for fname in fieldnames(T)
        fname = Meta.quot(fname)
        push!(cur.args, :(name === $fname))
        push!(cur.args, :(return $_literal_getproperty(sym, Val{$fname}())))
        push!(cur.args, Expr(:elseif))
        cur = cur.args[end]
    end
    cur.head = :block
    push!(cur.args, quote
        if @isdefined(FieldError)
            throw(FieldError($T, name))
        else
            error("type $($T) has no field $(name). Available fields are $($(fieldnames(T)))")
        end
    end)
    return chain
end

"""
    $TYPEDEF

Struct used as operation for symbolic getproperty on `SymStruct{T}` with field `field`.
"""
struct SymbolicGetproperty{T, field} end

"""
    field_name(f::SymbolicGetproperty{T, field}) -> Symbol

Return the name of the struct field accessed by `f`.
"""
field_name(::SymbolicGetproperty{T, field}) where {T, field} = field

function (f::SymbolicGetproperty{T})(x::SymbolicT) where {T}
    unwrap(f(SymStruct{T}(x)))
end
function (::SymbolicGetproperty{T, field})(x::SymStruct{T}) where {T, field}
    _literal_getproperty(x, Val{field}())
end
function (::SymbolicGetproperty{T, field})(x::T) where {T, field}
    getproperty(x, field)
end

function SymbolicUtils.promote_symtype(::SymbolicGetproperty{T, field}, x::SymbolicUtils.TypeT) where {T, field}
    @assert x == T
    fieldtype(x, field)
end

function SymbolicUtils.promote_shape(::SymbolicGetproperty{T, field},
                                     @nospecialize(x::SymbolicUtils.ShapeT)) where {T, field}
    @assert x isa SymbolicUtils.ShapeVecT && isempty(x)
    field_shape(T, Val{field}())
end

"""
    $TYPEDSIGNATURES

Called by the generated `getproperty` for `SymStruct`. Performs symbolic field access.
"""
function _literal_getproperty(sym::SymStruct{T}, ::Val{name}) where {T, name}
    fT = fieldtype(T, name)
    fShape = field_shape(T, Val{name}())
    fname = BSImpl.Const{VartypeT}(name)
    _struct = unwrap(sym)
    args = ArgsT{VartypeT}((_struct,))
    val = BSImpl.Term{VartypeT}(SymbolicGetproperty{T, name}(), args; type = fT, shape = fShape)
    if has_symwrapper(fT)
        return wrapper_type(fT)(val)
    else
        return val
    end
end

"""
    $TYPEDSIGNATURES

Obtain the shape of the value obtained by accessing field `name` of type `T`. Only
implemented by `@symstruct` via the `shape` option.
"""
function field_shape(::Type{T}, ::Val{name}) where {T, name}
    shape_from_type(fieldtype(T, name))
end

"""
    shape_from_type(::Type{T})

Infer the symbolic shape of a value of type `T`. For `AbstractArray{E, N}` subtypes,
returns `SymbolicUtils.Unknown(N)` (unknown shape with known rank). For all other types,
returns an empty `ShapeVecT`, indicating a scalar.
"""
shape_from_type(::Type{A}) where {T, N, A <: AbstractArray{T, N}} = SymbolicUtils.Unknown(N)
shape_from_type(::Type{T}) where {T} = SymbolicUtils.ShapeVecT()

function SymbolicUtils.show_call(io::IO, @nospecialize(f::SymbolicGetproperty), x::SymbolicT)
    fname = field_name(f)::Symbol
    @match x begin
        BSImpl.Term(; args) => print(io, args[1])
    end
    print(io, ".")
    print(io, fname)
end

function Base.show(io::IO, x::SymStruct)
    show(io, unwrap(x))
end

function SymbolicUtils.Code.function_to_expr(@nospecialize(f::SymbolicGetproperty), x::SymbolicT, st)
    out = get(st.rewrites, x, nothing)
    out === nothing  || return out

    fname = field_name(f)::Symbol
    args = @match x begin
        BSImpl.Term(; args) => args
    end
    return :($(SymbolicUtils.Code.toexpr(args[1], st)).$fname)
end

"""
    symstruct_field_supports_linear_indexing(::Type{T}, ::Val{field}) -> Bool

Return `true` if field `field` of struct type `T` can participate in linear indexing of a
`SymStruct{T}`. A field is indexable when its shape is fully known (i.e. not
`SymbolicUtils.Unknown`). Nested `SymStruct` fields are recursively checked via
[`symstruct_supports_linear_indexing`](@ref).
"""
function symstruct_field_supports_linear_indexing(::Type{T}, ::Val{field}) where {T, field}
    fT = fieldtype(T, field)
    if has_symwrapper(fT) && wrapper_type(fT) === SymStruct{fT}
        return symstruct_supports_linear_indexing(fT)
    end
    if fT <: Union{AbstractArray, Tuple}
        efT = eltype(fT)
        return field_shape(T, Val{field}()) isa SU.ShapeVecT && (
            !has_symwrapper(efT) || wrapper_type(efT) !== SymStruct{efT} ||
                symstruct_supports_linear_indexing(efT)
        )
    end
    return field_shape(T, Val{field}()) isa SU.ShapeVecT
end

"""
    symstruct_supports_linear_indexing(::Type{T}) -> Bool

Return `true` if every field of struct type `T` supports linear indexing, meaning that a
`SymStruct{T}` can be indexed with a single integer via `getindex`. This requires each
field (and any nested `SymStruct` fields) to have a fully known, statically determinable
shape.
"""
@generated function symstruct_supports_linear_indexing(::Type{T}) where {T}
    expr = Expr(:&&)
    for i in 1:fieldcount(T)
        fnameval = Val{fieldname(T, i)}()
        push!(expr.args, :($symstruct_field_supports_linear_indexing(T, $fnameval)))
    end
    return expr
end

"""
    symstruct_field_length(::Type{T}, ::Val{field}) -> Int

Return the number of scalar symbolic elements contributed by field `field` of struct type
`T` to the linear index space of a `SymStruct{T}`. Scalars contribute `1`, arrays
contribute the product of their shape ranges, and nested `SymStruct` fields contribute
their own [`symstruct_length`](@ref).
"""
function symstruct_field_length(::Type{T}, ::Val{field}) where {T, field}
    fT = fieldtype(T, field)
    if has_symwrapper(fT) && wrapper_type(fT) === SymStruct{fT}
        return symstruct_length(fT)
    end
    # The entire function is well-inferred, and this edge case
    # allows const-prop to collapse `symstruct_length` to a constant
    # when all fields are scalars.
    if !(fT <: Union{AbstractArray, Tuple})
        return 1
    end
    efT = eltype(fT)
    baselen = prod(length, field_shape(T, Val{field}())::SU.ShapeVecT; init = 1)
    if has_symwrapper(efT) && wrapper_type(efT) === SymStruct{efT}
        return baselen * symstruct_length(efT)
    end
    return baselen
end

"""
    symstruct_length(::Type{T}) -> Int

Return the total number of scalar symbolic elements in the linear index space of a
`SymStruct{T}`, equal to the sum of [`symstruct_field_length`](@ref) over all fields.
This is the value returned by `length` on a `SymStruct{T}` instance.
"""
@generated function symstruct_length(::Type{T}) where {T}
    expr = Expr(:call, (+))
    for i in 1:fieldcount(T)
        push!(expr.args, :($symstruct_field_length(T, Val{$(QuoteNode(fieldname(T, i)))}())))
    end
    return expr
end

Base.length(s::SymStruct{T}) where {T} = symstruct_length(T)

"""
    symstruct_checkbounds([Bool], s::SymStruct{T}, idx::Integer)

Check whether linear index `idx` is within bounds for `s`. The two-argument form (without
`Bool`) throws a `BoundsError` on failure. The three-argument form with `Bool` as first
argument returns `true` or `false` without throwing.
"""
function symstruct_checkbounds(::Type{Bool}, s::SymStruct{T}, idx::Integer) where {T}
    return 1 <= idx <= symstruct_length(T)
end
function symstruct_checkbounds(s::SymStruct{T}, idx::Integer) where {T}
    symstruct_checkbounds(Bool, s, idx) || throw(BoundsError(s, idx))
end

"""
    symstruct_field_range(::Type{T}, ::Val{field}) -> UnitRange{Int}

Return the range of linear indices within a `SymStruct{T}` that correspond to field
`field`. Always starts at `1`; the upper bound equals
[`symstruct_field_length`](@ref)`(T, Val{field}())`.
"""
function symstruct_field_range(::Type{T}, ::Val{field}) where {T, field}
    return 1:symstruct_field_length(T, Val{field}())
    fT = fieldtype(T, field)
    if has_symwrapper(fT) && wrapper_type(fT) === SymStruct{fT}
        return 1:symstruct_length(fT)
    end
    # Same fast-path as in `symstruct_field_length`
    if !(fT <: Union{AbstractArray, Tuple})
        return 1:1
    end
    fshape = field_shape(T, Val{field}())::SU.ShapeVecT
    return 1:prod(length, fshape; init = 1)
end

"""
    symstruct_field_getindex(s::SymStruct{T}, ::Val{field}, i::Integer) -> SymbolicT

Return the `i`-th scalar symbolic element within field `field` of `s`. For nested
`SymStruct` fields, delegates recursively. For array fields whose elements are nested
`SymStruct`s, maps `i` to the appropriate element and inner index.
"""
function symstruct_field_getindex(s::SymStruct{T}, ::Val{field}, i::Integer) where {T, field}
    fT = fieldtype(T, field)
    if has_symwrapper(fT) && wrapper_type(fT) === SymStruct{fT}
        return getindex(SymbolicGetproperty{T, field}()(s), i)
    end
    fval = SymbolicGetproperty{T, field}()(s)
    if fval isa Num
        return unwrap(fval)
    elseif fval isa SymbolicT
        return fval
    end
    fval = unwrap(fval)
    efT = eltype(fT)
    if has_symwrapper(efT) && wrapper_type(efT) === SymStruct{efT}
        N = symstruct_length(efT)
        return SymStruct{efT}(fval[SymbolicUtils.stable_eachindex(fval)[div(i - 1, N) + 1]])[mod1(i, N)]
    end
    return fval[SymbolicUtils.stable_eachindex(fval)[i]]
end

@generated function Base.getindex(s::SymStruct{T}, idx::Integer) where {T}
    conditions = Expr[]
    bodies = Expr[]
    offset = 0
    body = Expr(:block)
    for i in 1:fieldcount(T)
        fname = fieldname(T, i)
        ftype = fieldtype(T, i)
        fieldval = Val{fieldname(T, i)}()
        push!(body.args, :(if idx in $symstruct_field_range(T, $fieldval)
            return $symstruct_field_getindex(s, $fieldval, idx)::$SymbolicT
        else
            idx -= $symstruct_field_length(T, $fieldval)
        end))
    end

    return Expr(:block, Expr(:call, symstruct_checkbounds, :s, :idx), body, :(error("Unreachable")))
end

function Base.iterate(s::SymStruct{T}) where {T}
    return (s[1], 2)
end
function Base.iterate(s::SymStruct{T}, i::Int) where {T}
    i > length(s) && return nothing
    return (s[i], i + 1)
end
Base.eltype(s::SymStruct) = SymbolicT

