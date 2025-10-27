using SymbolicUtils: FnType, Sym, metadata
using Setfield

const IndexMap = Dict{Char,Char}(
            '-' => '₋',
            '0' => '₀',
            '1' => '₁',
            '2' => '₂',
            '3' => '₃',
            '4' => '₄',
            '5' => '₅',
            '6' => '₆',
            '7' => '₇',
            '8' => '₈',
            '9' => '₉')

abstract type AbstractVariableMetadata end
"""
    $TYPEDEF

Symbolic metadata key for storing the default value of a symbolic variable.
"""
struct VariableDefaultValue <: AbstractVariableMetadata end
"""
    $TYPEDEF

Symbolic metadata key for storing the macro used to create a symbolic variable.
"""
struct VariableSource <: AbstractVariableMetadata end

function setdefaultval(x, val)
    val === nothing && return x
    sh = shape(x)
    if sh isa SymbolicUtils.Unknown
        @assert sh.ndims == -1 || ndims(val) == sh.ndims """
        Variable $x must have default of matching `ndims`. Got $val with `ndims` \
        $(ndims(val)).
        """
    else
        @assert isempty(sh) || symtype(x) <: FnType || size(x) == size(val) """
        Variable $x must have default of matching size. Got $val with size \
        $(size(val)).
        """
    end
    setmetadata(x, VariableDefaultValue, val)
end

function map_subscripts(indices)
    str = string(indices)
    join(IndexMap[c] for c in str)
end


# Build variables more easily
"""
    $(TYPEDSIGNATURES)

Parse variables using the syntax expected by `@variables`. Used for implementing custom
macros similar to `@variables`. `macroname` refers to the name of the macro creating the
variables. This is stored in the `VariableSource` metadata of created variables. `type`
is the default type of created variables. `x` is the tuple of expressions passed to the
macro. `transform` is an optional function that takes constructed variables and performs
custom postprocessing to them, returning the created variables. This function returns the
`Expr` for constructing the parsed variables.
"""
function parse_vars(macroname, type, x, transform = identity)
    ex = Expr(:block)
    var_names = Expr(:vect)
    # if parsing things in the form of
    # begin
    #     x
    #     y
    #     z
    # end
    x = x isa Tuple && first(x) isa Expr && first(x).head == :tuple ? first(x).args : x # tuple handling
    x = flatten_expr!(x)
    cursor = 0
    isoption(ex) = Meta.isexpr(ex, [:vect, :vcat, :hcat])
    while cursor < length(x)
        cursor += 1
        var_expr = x[cursor]

        default = nothing
        options = nothing
        if Meta.isexpr(var_expr, :(=))
            var_expr, default = var_expr.args
            # defaults with metadata for function variables
            if Meta.isexpr(default, :block)
                Base.remove_linenums!(default)
                default = only(default.args)
            end
            if Meta.isexpr(default, :tuple) && length(default.args) == 2 && isoption(default.args[2])
                options = default.args[2].args
                default = default.args[1]
            end
            default = esc(default)
        end
        parse_result = SymbolicUtils.parse_variable(var_expr; default_type = type)
        handle_nonconcrete_symtype!(parse_result)
        sym = SymbolicUtils.sym_from_parse_result(parse_result, VartypeT)
        sym = handle_maybe_dependent_variable!(parse_result, sym, type)

        if options === nothing && cursor < length(x) && isoption(x[cursor + 1])
            options = x[cursor + 1].args
            cursor += 1
        end
        sym = _add_metadata(parse_result, sym, default, macroname, options)
        sym = handle_maybe_callandwrap!(parse_result, sym)
        sym = Expr(:call, transform, Expr(:call, wrap, sym))

        if parse_result[:isruntime]
            varname = Symbol(parse_result[:name])
        else
            varname = esc(parse_result[:name])
        end
        push!(var_names.args, varname)
        push!(ex.args, Expr(:(=), varname, sym))
    end
    push!(ex.args, var_names)
    return ex
end

function handle_nonconcrete_symtype!(parse_result)
    type = parse_result[:type]
    if type == :Complex
        parse_result[:type] = :(Complex{Real})
    end
    if Meta.isexpr(type, :curly)
        if type.args[1] in (:Array, :Vector, :Matrix, Array, Vector, Matrix) && type.args[2] == :Complex
            type.args[2] = :(Complex{Real})
        end
        if type.args[1] == :FnType || type.args[1] == SymbolicUtils.FnType
            if Meta.isexpr(type.args[2], :curly) # Tuple{...}
                for i in 2:length(type.args[2].args)
                    if type.args[2].args[i] == :Complex
                        type.args[2].args[i] = :(Complex{Real})
                    end
                end
            end
            if type.args[3] == :Complex
                type.args[3] = :(Complex{Real})
            end
            for parse_arg in parse_result[:args]
                handle_nonconcrete_symtype!(parse_arg)
            end
        end
    end
    return nothing
end

function parse_result_is_dependent_variable(parse_result)
    # This means it is a function call
    return haskey(parse_result, :args) &&
    # This checks `fnT` in `FnType{argsT, retT, fnT}`
        parse_result[:type].args[4] === Nothing &&
    # This ensures all arguments have defined names
        all(n -> n !== nothing && n != :..,
            (get(arg, :name, nothing) for arg in parse_result[:args]))
end

function handle_maybe_dependent_variable!(parse_result, sym, type)
    # is a function call and the function doesn't have a type and all arguments
    # are named
    parse_result_is_dependent_variable(parse_result) || return sym

    args = parse_result[:args]
    argnames = Any[get(arg, :name, nothing) for arg in args]
    # if the last arg is a `Vararg`, splat it
    if !isempty(args) && Meta.isexpr(args[end][:type], :curly) && args[end][:type].args[1] == :Vararg
        argnames[end] = Expr(:..., argnames[end])
    end
    # Turn the result into something of the form `@variables x(..)`.
    # This makes it so that the `FnType` is recognized as a dependent variable
    # according to `SymbolicUtils.is_function_symtype`
    parse_result[:args] = [SymbolicUtils.parse_variable(:(..); default_type = type)]
    parse_result[:type].args[2] = Tuple
    # Re-create the `Sym`
    sym = SymbolicUtils.sym_from_parse_result(parse_result, VartypeT)
    # Change the type
    parse_result[:type] = parse_result[:type].args[3]
    # Call the `Sym` with the arguments to create a dependent variable.
    map!(esc, argnames, argnames)
    sym = Expr(:call, sym)
    append!(sym.args, argnames)
    return sym
end

function handle_maybe_callandwrap!(parse_result, sym)
    type = parse_result[:type]
    if Meta.isexpr(type, :curly) && (type.args[1] == :FnType || type.args[1] === SymbolicUtils.FnType)
        sym = Expr(:call, CallAndWrap, sym)
    end
    return sym
end

function _add_metadata(parse_result, var::Expr, default, macroname::Symbol, metadata::Union{Nothing, Vector{Any}})
    @nospecialize var default metadata
    if default !== nothing
        var = Expr(:call, setdefaultval, var, default)
    end
    varname = parse_result[:name]
    if parse_result[:isruntime]
        varname = esc(varname)
    else
        varname = Meta.quot(varname)
    end
    var = Expr(:call, setmetadata, var, VariableSource, Expr(:tuple, Meta.quot(macroname), varname))
    metadata === nothing && return var
    for ex in metadata
        Meta.isexpr(ex, :(=)) || error("Metadata must of the form of `key = value`")
        key, value = ex.args
        key_type = option_to_metadata_type(Val{key}())::DataType
        var = Expr(:call, setmetadata, var, key_type, esc(value))
    end
    return var
end

"""
    $(TYPEDSIGNATURES)

Define a new metadata key assignable in `@variables`. This function should take `Val{name}`
where `name` is a `Symbol`, and return the key type for the given metadata name `name`. For
example,

```julia
Symbolics.option_to_metadata_type(::Val{:custom_name}) = CustomType
```

Allows the following syntax:

```julia
@variables x [custom_name = 1]
```

And stores `1` as the value associated with the `CustomType` key in the symbolic metadata
of `x`.
"""
function option_to_metadata_type(::Val{opt}) where {opt}
    throw(Base.Meta.ParseError("unknown property type $opt"))
end

# add enough additional methods that the compiler gives up on specializing this
# and downstream definitions don't cause massive invalidation.
option_to_metadata_type(::Val{:_____!_internal_1}) = error("Invalid option")
option_to_metadata_type(::Val{:_____!_internal_2}) = error("Invalid option")
option_to_metadata_type(::Val{:_____!_internal_3}) = error("Invalid option")
option_to_metadata_type(::Val{:_____!_internal_4}) = error("Invalid option")
option_to_metadata_type(::Val{:_____!_internal_5}) = error("Invalid option")

"""

Define one or more unknown variables.

```julia
@variables t α σ(..) β[1:2]
@variables w(..) x(t) y z(t, α, x)

expr = β[1]* x + y^α + σ(3) * (z - t) - β[2] * w(t - 1)
```

`(..)` signifies that the value should be left uncalled.

Symbolics supports creating variables that denote an array of some size.

```jldoctest
julia> @variables x[1:3]
1-element Vector{Symbolics.Arr{Num, 1}}:
 x[1:3]

julia> @variables y[1:3, 1:6] # support for  tensors
1-element Vector{Symbolics.Arr{Num, 2}}:
 y[1:3,1:6]

julia> @variables t z(t)[1:3] # also works for dependent variables
2-element Vector{Any}:
 t
  (z(t))[1:3]
```

A symbol or expression that represents an array can be turned into an array of
symbols or expressions using the `scalarize` function.

```jldoctest
julia> Symbolics.scalarize(z)
3-element Vector{Num}:
 (z(t))[1]
 (z(t))[2]
 (z(t))[3]
```

Note that `@variables` returns a vector of all the defined variables.

`@variables` can also take runtime symbol values by the `\$` interpolation
operator, and in this case, `@variables` doesn't automatically assign the value,
instead, it only returns a vector of symbolic variables. All the rest of the
syntax also applies here.

```jldoctest
julia> a, b, c = :runtime_symbol_value, :value_b, :value_c
(:runtime_symbol_value, :value_b, :value_c)

julia> vars = @variables t \$a \$b(t) \$c(t)[1:3]
4-element Vector{Any}:
      t
 runtime_symbol_value
   value_b(t)
       (value_c(t))[1:3]

julia> (t, a, b, c)
(t, :runtime_symbol_value, :value_b, :value_c)
```
"""
macro variables(xs...)
    parse_vars(:variables, Real, xs)
end

const _fail = Dict()

getsource(x, val=_fail) = getmetadata(unwrap(x), VariableSource, val)

SymbolicIndexingInterface.symbolic_type(::Type{Symbolics.Num}) = ScalarSymbolic()
SymbolicIndexingInterface.symbolic_type(::Type{Symbolics.Arr{T, N}}) where {T, N} = ArraySymbolic()

SymbolicIndexingInterface.hasname(x::Union{Num,Arr,Complex{Num}}) = hasname(unwrap(x))
function SymbolicIndexingInterface.getname(x::Union{Num, Arr, Complex{Num}})
    SymbolicIndexingInterface.getname(unwrap(x))
end

function SymbolicIndexingInterface.symbolic_evaluate(ex::Union{Num, Arr, BasicSymbolic, Equation, Inequality}, d::Dict; kwargs...)
    val = fixpoint_sub(ex, d; fold = Val(true), kwargs...)
    return _recursive_unwrap(val)
end

for T in [LinearAlgebra.UpperTriangular, LinearAlgebra.LowerTriangular]
    @eval function _recursive_unwrap(val::$T)
        $T(_recursive_unwrap(collect(val)))
    end
end

function _recursive_unwrap(val)
    if symbolic_type(val) == NotSymbolic() && val isa Union{AbstractArray, Tuple}
        if parent(val) !== val
            return Setfield.@set val.parent = _recursive_unwrap(parent(val))
        end
        return _recursive_unwrap.(val)
    else
        return unwrap(val)
    end
end

function _recursive_unwrap(val::AbstractSparseArray)
    if val isa AbstractSparseVector
        (Is, Vs) = findnz(val)
        Vs = _recursive_unwrap.(Vs)
        return SparseVector(length(val), Is, Vs)
    else
        (Is, Js, Vs) = findnz(val)
        Vs = _recursive_unwrap.(Vs)
        return sparse(Is, Js, Vs, size(val)...) 
    end
end

struct FPSubFilterer{O} end

function (::FPSubFilterer{O})(ex::BasicSymbolic{T}) where {T, O}
    @match ex begin
        BSImpl.Term(; f) && if f isa Operator end => !(f isa O)
        _ => SymbolicUtils.default_substitute_filter(ex)
    end
end

"""
    fixpoint_sub(expr, dict; operator = Nothing, maxiters = 1000)

Given a symbolic expression, equation or inequality `expr` perform the substitutions in
`dict` recursively until the expression does not change. Substitutions that depend on one
another will thus be recursively expanded. For example,
`fixpoint_sub(x, Dict(x => y, y => 3))` will return `3`. The `operator` keyword can be
specified to prevent substitution of expressions inside operators of the given type. The
`maxiters` keyword is used to limit the number of times the substitution can occur to avoid
infinite loops in cases where the substitutions in `dict` are circular
(e.g. `[x => y, y => x]`).
"""
function fixpoint_sub(x, dict; operator = Nothing, maxiters = 1000, kw...)
    y = substitute(x, dict; filterer=FPSubFilterer{operator}(), kw...)
    iters = maxiters
    while !isequal(x, y) && iters > 0
        y = x
        x = substitute(y, dict; filterer=FPSubFilterer{operator}(), kw...)
        iters -= 1
    end

    if !isequal(x, y)
        @warn "Did not converge after `maxiters = $maxiters` substitutions. Either there is a cycle in the rules or `maxiters` needs to be higher."
    end

    return x
end
function fixpoint_sub(x::SparseMatrixCSC, dict; operator = Nothing, maxiters = 1000, kw...)
    I, J, V = findnz(x)
    V = fixpoint_sub(V, dict; operator, maxiters, kw...)
    m, n = size(x)
    return sparse(I, J, V, m, n)
end

function is_array_of_symbolics(x)
    symbolic_type(x) == ArraySymbolic() && return true
    symbolic_type(x) == ScalarSymbolic() && return false
    x isa AbstractArray &&
        any(y -> symbolic_type(y) != NotSymbolic() || is_array_of_symbolics(y), x)
end

function getdefaultval(x, val=_fail)
    x = unwrap(x)
    val = getmetadata(x, VariableDefaultValue, val)
    if val !== _fail
        return val
    else
        @match x begin
            BSImpl.Term(; f, args) && if f === getindex end => begin
                idxs = Iterators.map(unwrap_const, Iterators.drop(args, 1))
                return getdefaultval(args[1], val)[idxs...]
            end
            _ => error("$x has no default value")
        end
    end
end

"""
    variables(name::Symbol, indices...)

Create a multi-dimensional array of individual variables named with subscript
notation. Use `@variables` instead to create symbolic array variables (as
opposed to array of variables). See `variable` to create one variable with
subscripts.

```jldoctest
julia> Symbolics.variables(:x, 1:3, 3:6)
3×4 Matrix{Num}:
 x₁ˏ₃  x₁ˏ₄  x₁ˏ₅  x₁ˏ₆
 x₂ˏ₃  x₂ˏ₄  x₂ˏ₅  x₂ˏ₆
 x₃ˏ₃  x₃ˏ₄  x₃ˏ₅  x₃ˏ₆
```
"""
function variables(name, indices...; T=Real)
    [variable(name, ij...; T=T) for ij in Iterators.product(indices...)]
end

"""
    variable(name::Symbol, idx::Integer...; T=Real)

Create a variable with the given name along with subscripted indices with the
`symtype=T`. When `T=FnType`, it creates a symbolic function.

```jldoctest
julia> Symbolics.variable(:x, 4, 2, 0)
x₄ˏ₂ˏ₀

julia> Symbolics.variable(:x, 4, 2, 0, T=Symbolics.FnType)
x₄ˏ₂ˏ₀⋆
```

Also see `variables`.
"""
function variable(name, idx...; T=Real)
    name_ij = Symbol(name, join(map_subscripts.(idx), "ˏ"))
    v = Sym{VartypeT}(name_ij; type = T)
    wrap(setmetadata(v, VariableSource, (:variables, name_ij)))
end

##### Renaming #####
# getname
# rename
# getindex parent
# calls
# symbolic function x[1:3](..)
#
# x_t
# sys.x

function renamed_metadata(metadata::Union{Nothing, SymbolicUtils.MetadataT}, name::Symbol)
    @nospecialize metadata
    if metadata === nothing
        return metadata
    elseif metadata isa Base.ImmutableDict{DataType, Any}
        newmeta = Base.ImmutableDict{DataType, Any}()
        for (k, v) in metadata
            if k === VariableSource
                v = v::NTuple{2, Symbol}
                v = (v[1], name)
            end
            newmeta = Base.ImmutableDict(newmeta, k, v)
        end
        return newmeta
    end
    error()
end

rename(x::Union{Num, Arr}, name) = wrap(rename(unwrap(x), name))

function rename(x::BasicSymbolic{T}, newname::Symbol) where {T}
    @match x begin
        BSImpl.Sym(; name, type, shape, metadata) => begin
            metadata = renamed_metadata(metadata, newname)
            return BSImpl.Sym{T}(newname; type, shape, metadata)
        end
        BSImpl.Term(; f, args, type, shape, metadata) && if f === getindex end => begin
            newargs = copy(parent(args))
            newargs[1] = rename(newargs[1], newname)
            return BSImpl.Term{T}(f, newargs; type, shape, metadata)
        end
        BSImpl.Term(; f, args, type, shape, metadata) && if f isa BasicSymbolic{T} end => begin
            f = rename(f, newname)
            metadata = renamed_metadata(metadata, newname)
            return BSImpl.Term{T}(f, args; type, shape, metadata)
        end
        _ => error("Cannot rename $x.")
    end
end

#### Callable ###

struct CallAndWrap{T}
    f::BasicSymbolic{VartypeT}
end

rettype_from_fntype(::Type{T}) where {A, R, T <: SymbolicUtils.FnType{A, R}} = R

function CallAndWrap(f::BasicSymbolic{VartypeT})
    @assert symtype(f) <: SymbolicUtils.FnType
    R = rettype_from_fntype(symtype(f))
    if hasmethod(wrapper_type, Tuple{Type{R}})
        CallAndWrap{wrapper_type(R)}(f)
    else
        f
    end
end

function (caw::CallAndWrap{T})(args...) where {T}
    T(caw.f(args...))
end

SymbolicIndexingInterface.symbolic_type(::Type{CallAndWrap{T}}) where {T} = ScalarSymbolic()

Base.isequal(a::CallAndWrap, b::CallAndWrap) = isequal(a.f,  b.f)
Base.isequal(a::BasicSymbolic{VartypeT}, b::CallAndWrap) = isequal(a,  b.f)
Base.isequal(a::CallAndWrap, b::BasicSymbolic{VartypeT}) = isequal(a.f,  b)
Base.hash(x::CallAndWrap, h::UInt) = hash(x.f, h)

function Base.show(io::IO, caw::CallAndWrap)
    show(io, caw.f)
    print(io, "⋆")
end

has_symwrapper(::Type{T}) where {A, R, T <: SymbolicUtils.FnType{A, R}} = has_symwrapper(R)
wrapper_type(::Type{T}) where {A, R, T <: SymbolicUtils.FnType{A, R}} = CallAndWrap{wrapper_type(R)}
is_wrapper_type(::Type{T}) where {T <: CallAndWrap} = true
wraps_type(::Type{T}) where {W, T <: CallAndWrap{W}} = FnType{A, R} where {A, R <: wraps_type(W)}
iswrapped(::CallAndWrap) = true
SymbolicUtils.unwrap(x::CallAndWrap) = x.f
SymbolicUtils.symtype(x::CallAndWrap) = symtype(x.f)
SymbolicIndexingInterface.getname(x::CallAndWrap) = getname(x.f)
SymbolicIndexingInterface.hasname(x::CallAndWrap) = hasname(x.f)
Symbolics.rename(x::CallAndWrap, name) = CallAndWrap(rename(x.f, name))
SymbolicUtils.getmetadata(x::CallAndWrap, args...) = SymbolicUtils.getmetadata(x.f, args...)
SymbolicUtils.setmetadata(x::CallAndWrap, args...) = SymbolicUtils.setmetadata(x.f, args...)
SymbolicUtils.hasmetadata(x::CallAndWrap, args...) = SymbolicUtils.hasmetadata(x.f, args...)
