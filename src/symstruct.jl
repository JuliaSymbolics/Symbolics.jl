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

SymbolicUtils.unwrap(x::SymStruct) = getfield(x, 1)
SymbolicUtils.infer_vartype(::Type{SymStruct{T}}) where {T} = VartypeT

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

function SymbolicUtils.promote_type(::SymbolicGetproperty{T, field}, x::SymbolicUtils.TypeT) where {T, field}
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
    args = ArgsT{VartypeT}((_struct, fname))
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
