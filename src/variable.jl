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

struct VariableDefaultValue end
struct VariableSource end

function recurse_and_apply(f, x)
    if symtype(x) <: AbstractArray
        getindex_posthook(x) do r,x,i...
            recurse_and_apply(f, r)
        end
    else
        f(x)
    end
end

function setdefaultval(x, val)
    if symtype(x) <: AbstractArray
        if val isa AbstractArray
            getindex_posthook(x) do r,x,i...
                setdefaultval(r, val[i...])
            end
        else
            getindex_posthook(x) do r,x,i...
                setdefaultval(r, val)
            end
        end
    else
        setmetadata(x,
                    VariableDefaultValue,
                    val)
    end
end

struct GetindexParent end

function scalarize_getindex(x, parent=x)
    if symtype(x) <: AbstractArray
        getindex_posthook(x) do r,x,i...
            scalarize_getindex(r, parent)
        end
    else
        xx = scalarize(x)
        xx = metadata(xx, metadata(x))
        if symtype(xx) <: FnType
            setmetadata(CallWithMetadata(xx, metadata(xx)), GetindexParent, parent)
        else
            setmetadata(xx, GetindexParent, parent)
        end
    end
end


function map_subscripts(indices)
    str = string(indices)
    join(IndexMap[c] for c in str)
end


function unwrap_runtime_var(v)
    isruntime = Meta.isexpr(v, :$) && length(v.args) == 1
    isruntime && (v = v.args[1])
    return isruntime, v
end

# Build variables more easily
function _parse_vars(macroname, type, x, transform=identity)
    ex = Expr(:block)
    var_names = Symbol[]
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
        v = x[cursor]

        # We need lookahead to the next `v` to parse
        # `@variables x [connect=Flow,unit=u]`
        nv = cursor < length(x) ? x[cursor+1] : nothing
        val = unit = connect = options = nothing

        # x = 1, [connect = flow; unit = u"m^3/s"]
        if Meta.isexpr(v, :(=))
            v, val = v.args
            if Meta.isexpr(val, :tuple) && length(val.args) == 2 && isoption(val.args[2])
                options = val.args[2].args
                val = val.args[1]
            end
        end

        type′ = type

        if Meta.isexpr(v, :(::))
            v, type′ = v.args
            type′ = type′ === :Complex ? Complex{type} : type′
        end


        # x [connect = flow; unit = u"m^3/s"]
        if isoption(nv)
            options = nv.args
            cursor += 1
        end

        isruntime, v = unwrap_runtime_var(v)
        iscall = Meta.isexpr(v, :call)
        isarray = Meta.isexpr(v, :ref)
        issym  = v isa Symbol
        @assert iscall || isarray || issym "@$macroname expects a tuple of expressions or an expression of a tuple (`@$macroname x y z(t) v[1:3] w[1:2,1:4]` or `@$macroname x y z(t) v[1:3] w[1:2,1:4] k=1.0`)"

        if iscall
            isruntime, fname = unwrap_runtime_var(v.args[1])
            call_args = map(last∘unwrap_runtime_var, @view v.args[2:end])
            var_name, expr = construct_vars(macroname, fname, type′, call_args, val, options, transform, isruntime)
        else
            var_name, expr = construct_vars(macroname, v, type′, nothing, val, options, transform, isruntime)
        end

        push!(var_names, var_name)
        push!(ex.args, expr)
    end
    rhs = build_expr(:vect, var_names)
    push!(ex.args, rhs)
    return ex
end

function construct_vars(macroname, v, type, call_args, val, prop, transform, isruntime)
    issym  = v isa Symbol
    isarray = isa(v, Expr) && v.head == :ref
    if isarray
        var_name = v.args[1]
        if Meta.isexpr(var_name, :(::))
            var_name, type′ = var_name.args
            type = type′ === :Complex ? Complex{type} : type′
        end
        isruntime, var_name = unwrap_runtime_var(var_name)
        indices = v.args[2:end]
        expr = _construct_array_vars(macroname, isruntime ? var_name : Meta.quot(var_name), type, call_args, val, prop, indices...)
    else
        var_name = v
        if Meta.isexpr(v, :(::))
            var_name, type′ = v.args
            type = type′ === :Complex ? Complex{type} : type′
        end
        expr = construct_var(macroname, isruntime ? var_name : Meta.quot(var_name), type, call_args, val, prop)
    end
    lhs = isruntime ? gensym(var_name) : var_name
    rhs = :($transform($expr))
    lhs, :($lhs = $rhs)
end

function option_to_metadata_type(::Val{opt}) where {opt}
    throw(Base.Meta.ParseError("unknown property type $opt"))
end

function setprops_expr(expr, props, macroname, varname)
    expr = :($setmetadata($expr, $VariableSource, ($(Meta.quot(macroname)), $varname,)))
    isnothing(props) && return expr
    for opt in props
        if !Meta.isexpr(opt, :(=))
            throw(Base.Meta.ParseError(
                "Variable properties must be in " *
                "the form of `a = b`. Got $opt."))
        end

        lhs, rhs = opt.args

        @assert lhs isa Symbol "the lhs of an option must be a symbol"
        expr = :($setmetadata($expr,
                              $(option_to_metadata_type(Val{lhs}())),
                       $rhs))
    end
    expr
end

struct CallWithMetadata{T,M} <: Symbolic{T}
    f::Symbolic{T}
    metadata::M
end

for f in [:istree, :operation, :arguments]
    @eval SymbolicUtils.$f(x::CallWithMetadata) = $f(x.f)
end

SymbolicUtils.Code.toexpr(x::CallWithMetadata, st) = SymbolicUtils.Code.toexpr(x.f, st)

CallWithMetadata(f) = CallWithMetadata(f, nothing)

function Base.show(io::IO, c::CallWithMetadata)
    show(io, c.f)
    print(io, "⋆")
end

function (f::CallWithMetadata)(args...)
    wrap(metadata(f.f(args...), metadata(f)))
end

function construct_var(macroname, var_name, type, call_args, val, prop)
    expr = if call_args === nothing
        :($Sym{$type}($var_name))
    elseif !isempty(call_args) && call_args[end] == :..
        :($CallWithMetadata($Sym{$FnType{Tuple, $type}}($var_name)))
    else
        :($Sym{$FnType{NTuple{$(length(call_args)), Any}, $type}}($var_name)($(map(x->:($value($x)), call_args)...)))
    end

    if val !== nothing
        expr = :($setdefaultval($expr, $val))
    end

    :($wrap($(setprops_expr(expr, prop, macroname, var_name))))
end

struct CallWith
    args
end

(c::CallWith)(f) = unwrap(f(c.args...))

function _construct_array_vars(macroname, var_name, type, call_args, val, prop, indices...)
    # TODO: just use Sym here
    ndim = length(indices)

    need_scalarize = false
    expr = if call_args === nothing
        ex = :($Sym{Array{$type, $ndim}}($var_name))
        :($setmetadata($ex, $ArrayShapeCtx, ($(indices...),)))
    elseif !isempty(call_args) && call_args[end] == :..
        need_scalarize = true
        ex = :($Sym{Array{$FnType{Tuple, $type}, $ndim}}($var_name))
        ex = :($setmetadata($ex, $ArrayShapeCtx, ($(indices...),)))
        :($map($CallWithMetadata, $ex))
    else
        # [(R -> R)(R) ....]
        need_scalarize = true
        ex = :($Sym{Array{$FnType{Tuple, $type}, $ndim}}($var_name))
        ex = :($setmetadata($ex, $ArrayShapeCtx, ($(indices...),)))
        :($map($CallWith(($(call_args...),)), $ex))
    end

    if val !== nothing
        expr = :($setdefaultval($expr, $val))
    end
    expr = setprops_expr(expr, prop, macroname, var_name)
    if need_scalarize
        expr = :($scalarize_getindex($expr))
    end

    expr = :($wrap($expr))

    return expr
end


"""

Define one or more unknown variables.

```julia
@variables t α σ(..) β[1:2]
@variables w(..) x(t) y z(t, α, x)

expr = β[1]* x + y^α + σ(3) * (z - t) - β[2] * w(t - 1)
```

`(..)` signifies that the value should be left uncalled.

Symbolics supports creating variables that denote an array of some size.

```julia
julia> @variables x[1:3]
1-element Vector{Symbolics.Arr{Num, 1}}:
 x[1:3]

julia> @variables y[1:3, 1:6] # support for  tensors
1-element Vector{Symbolics.Arr{Num, 2}}:
 y[1:3,1:6]

julia> @variables t z[1:3](t) # also works for dependent variables
2-element Vector{Any}:
 t
  (map(#5, z))[1:3]
```

A symbol or expression that represents an array can be turned into an array of
symbols or expressions using the `scalarize` function.

```julia
julia> Symbolics.scalarize(z)
3-element Vector{Num}:
 z[1](t)
 z[2](t)
 z[3](t)
```

Note that `@variables` returns a vector of all the defined variables.

`@variables` can also take runtime symbol values by the `\$` interpolation
operator, and in this case, `@variables` doesn't automatically assign the value,
instead, it only returns a vector of symbolic variables. All the rest of the
syntax also applies here.

```julia
julia> a, b, c = :runtime_symbol_value, :value_b, :value_c
:runtime_symbol_value

julia> vars = @variables t \$a \$b(t) \$c[1:3](t)
4-element Vector{Any}:
      t
 runtime_symbol_value
   value_b(t)
       (map(#9, value_c))[1:3]

julia> (t, a, b, c)
(t, :runtime_symbol_value, :value_b, :value_c)
```
"""
macro variables(xs...)
    esc(_parse_vars(:variables, Real, xs))
end

TreeViews.hastreeview(x::Sym) = true

function TreeViews.treelabel(io::IO,x::Sym,
                             mime::MIME"text/plain" = MIME"text/plain"())
  show(io,mime,Text(x.name))
end

struct Namespace{T} <: Symbolic{T}
    parent::Any
    named::Symbolic{T}
    function Namespace(p, n)
        n isa Namespace && error("Ill-formed namespacing. $n shouldn't be a namespace.")
        new{symtype(n)}(p, n)
    end
end

Base.hash(ns::Namespace, salt::UInt) = hash(ns.named, hash(getname(ns.parent), salt ⊻ 0x906e89687f904e4a))
SymbolicUtils.metadata(ns::Namespace) = SymbolicUtils.metadata(ns.named)
SymbolicUtils.metadata(ns::Namespace, meta) = @set ns.named = SymbolicUtils.metadata(ns.named, meta)
SymbolicUtils.setmetadata(ns::Namespace, typ::DataType, data) = @set ns.named = SymbolicUtils.setmetadata(ns.named, typ, data)
Base.nameof(x::Namespace) = getname(x)
function SymbolicUtils.Code.toexpr(ns::Namespace, st)
    if haskey(st.symbolify, ns)
        st.symbolify[ns]
    else
        :($(getname(ns.parent)).$(SymbolicUtils.Code.toexpr(ns.named, st)))
    end
end
Base.show(io::IO, x::Namespace) = print(io, getname(x))
function Base.isequal(x::Namespace, y::Namespace)
    isequal(x.named, y.named) && (x.parent === y.parent || isequal(x.parent, y.parent))
end
# Namespace must be treated as a variable
SymbolicUtils.Code.get_symbolify(ns::Namespace) = (ns,)

const _fail = Dict()

_getname(x, _) = nameof(x)
_getname(x::Symbol, _) = x
function _getname(x::Symbolic, val)
    if istree(x) && (op = operation(x)) isa Namespace
        return getname(op)
    end
    ss = getsource(x, nothing)
    if ss === nothing
        ss = getsource(getparent(x), val)
    end
    ss === _fail && throw(ArgumentError("Variable $x doesn't have a source defined."))
    ss[2]
end
_getname(x::Namespace, val) = Symbol(getname(x.parent), :(.), getname(x.named, val))

getsource(x, val=_fail) = getmetadata(unwrap(x), VariableSource, val)

getname(x, val=_fail) = _getname(unwrap(x), val)

function getparent(x, val=_fail)
    maybe_parent = getmetadata(x, Symbolics.GetindexParent, nothing)
    if maybe_parent !== nothing
        return maybe_parent
    else
        if istree(x) && operation(x) === getindex
            return arguments(x)[1]
        end
    end
    val === _fail && throw(ArgumentError("Cannot find the parent of $x."))
    return val
end

function getdefaultval(x, val=_fail)
    x = unwrap(x)
    val = getmetadata(x, VariableDefaultValue, val)
    if val !== _fail
        return val
    else
        error("$x has no default value")
    end
end

"""
    variables(name::Symbol, indices...)

Create a multi-dimensional array of individual variables named with subscript
notation. Use `@variables` instead to create symbolic array variables (as
opposed to array of variables). See `variable` to create one variable with
subscripts.

```julia-repl
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

```julia-repl
julia> Symbolics.variable(:x, 4, 2, 0)
x₄ˏ₂ˏ₀

julia> Symbolics.variable(:x, 4, 2, 0, T=Symbolics.FnType)
x₄ˏ₂ˏ₀⋆
```

Also see `variables`.
"""
function variable(name, idx...; T=Real)
    name_ij = Symbol(name, join(map_subscripts.(idx), "ˏ"))
    if T <: FnType
        first(@variables $name_ij(..))
    else
        first(@variables $name_ij::T)
    end
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

function rename_getindex_source(x, parent=x)
    getindex_posthook(x) do r,x,i...
        hasmetadata(r, GetindexParent) ? setmetadata(r, GetindexParent, parent) : r
    end
end

function rename_source(from, to)
    if hasmetadata(from, VariableSource)
        s = getmetadata(from, VariableSource)
        setmetadata(to, VariableSource, (s[1], getname(to)))
    else
        to
    end
end

function rename(x::Sym, name)
    xx = @set! x.name = name
    xx = rename_source(x, xx)
end

rename(x::Union{Num, Arr}, name) = wrap(rename(unwrap(x), name))
function rename(x::ArrayOp, name)
    t = x.term
    args = arguments(t)
    @show x
    # Hack:
    #@assert operation(t) === (map) && (args[1] isa CallWith || args[1] == CallWithMetadata)
    rn = rename(args[2], name)

    xx = metadata(operation(t)(args[1], rn),
                  metadata(x))
    rename_getindex_source(rename_source(x, xx))
end

function rename(x::CallWithMetadata, name)
    rename_source(x, CallWithMetadata(rename(x.f, name), x.metadata))
end

function rename(x::Symbolic, name)
    if operation(x) isa Sym
        @assert x isa Term
        @set! x.f = rename(operation(x), name)
        @set! x.hash = Ref{UInt}(0)
        return x
    else
        error("can't rename $x to $name")
    end
end

# Deprecation below

struct Variable{T} end

function (::Type{Variable{T}})(s, i...) where {T}
    Base.depwarn("Variable{T}(name, idx...) is deprecated, use variable(name, idx...; T=T)", :Variable, force=true)
    variable(s, i...; T=T)
end

(::Type{Variable})(s, i...) = Variable{Real}(s, i...)

function (::Type{Sym{T}})(s, x, i...) where {T}
    Base.depwarn("Sym{T}(name, x, idx...) is deprecated, use variable(name, x, idx...; T=T)", :Variable, force=true)
    variable(s, x, i...; T=T)
end
(::Type{Sym})(s, x, i...) = Sym{Real}(s, x, i...)
