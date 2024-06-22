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
struct VariableDefaultValue <: AbstractVariableMetadata end
struct VariableSource <: AbstractVariableMetadata end

function recurse_and_apply(f, x)
    if symtype(x) <: AbstractArray
        getindex_posthook(x) do r,x,i...
            recurse_and_apply(f, r)
        end
    else
        f(x)
    end
end

function set_scalar_metadata(x, V, val)
    if symtype(x) <: AbstractArray
        x = if val isa AbstractArray
            getindex_posthook(x) do r,x,i...
                set_scalar_metadata(r, V, val[i...])
            end
        else
            getindex_posthook(x) do r,x,i...
                set_scalar_metadata(r, V, val)
            end
        end
    end
    setmetadata(x, V, val)
end
setdefaultval(x, val) = set_scalar_metadata(x, VariableDefaultValue, val)

struct GetindexParent end

function scalarize_getindex(x, parent=Ref{Any}(x))
    if symtype(x) <: AbstractArray
        parent[] = getindex_posthook(x) do r,x,i...
            scalarize_getindex(r, parent)
        end
    else
        xx = unwrap(scalarize(x))
        xx = metadata(xx, metadata(x))
        if symtype(xx) <: FnType
            setmetadata(CallWithMetadata(xx, metadata(xx)), GetindexParent, parent[])
        else
            setmetadata(xx, GetindexParent, parent[])
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
        if iscall && Meta.isexpr(v.args[1], :ref)
            @warn("The variable syntax $v is deprecated. Use $(Expr(:ref, Expr(:call, v.args[1].args[1], v.args[2]), v.args[1].args[2:end]...)) instead.
                  The former creates an array of functions, while the latter creates an array valued function.
                  The deprecated syntax will cause an error in the next major release of Symbolics.
                  This change will facilitate better implementation of various features of Symbolics.")
        end
        issym  = v isa Symbol
        @assert iscall || isarray || issym "@$macroname expects a tuple of expressions or an expression of a tuple (`@$macroname x y z(t) v[1:3] w[1:2,1:4]` or `@$macroname x y z(t) v[1:3] w[1:2,1:4] k=1.0`)"

        if isarray && Meta.isexpr(v.args[1], :call)
            # This is the new syntax
            isruntime, fname = unwrap_runtime_var(v.args[1].args[1])
            call_args = map(last∘unwrap_runtime_var, @view v.args[1].args[2:end])
            size = v.args[2:end]
            var_name, expr = construct_dep_array_vars(macroname, fname, type′, call_args, size, val, options, transform, isruntime)
        elseif iscall
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

function construct_dep_array_vars(macroname, lhs, type, call_args, indices, val, prop, transform, isruntime)
    ndim = :($length(($(indices...),)))
    vname = !isruntime ? Meta.quot(lhs) : lhs
    if call_args[1] == :..
        ex = :($CallWithMetadata($Sym{$FnType{Tuple, Array{$type, $ndim}}}($vname)))
    else
        ex = :($Sym{$FnType{Tuple, Array{$type, $ndim}}}($vname)(map($unwrap, ($(call_args...),))...))
    end
    ex = :($setmetadata($ex, $ArrayShapeCtx, ($(indices...),)))

    if val !== nothing
        ex = :($setdefaultval($ex, $val))
    end
    ex = setprops_expr(ex, prop, macroname, Meta.quot(lhs))
    #ex = :($scalarize_getindex($ex))

    ex = :($wrap($ex))

    if call_args[1] == :..
        ex = :($transform($ex))
    end
    if isruntime
        lhs = gensym(lhs)
    end
    lhs, :($lhs = $ex)
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
        expr = :($set_scalar_metadata($expr,
                              $(option_to_metadata_type(Val{lhs}())),
                       $rhs))
    end
    expr
end

struct CallWithMetadata{T,M} <: Symbolic{T}
    f::Symbolic{T}
    metadata::M
end

for f in [:iscall, :operation, :arguments]
    @eval SymbolicUtils.$f(x::CallWithMetadata) = $f(x.f)
end

SymbolicUtils.Code.toexpr(x::CallWithMetadata, st) = SymbolicUtils.Code.toexpr(x.f, st)

CallWithMetadata(f) = CallWithMetadata(f, nothing)

function Base.show(io::IO, c::CallWithMetadata)
    show(io, c.f)
    print(io, "⋆")
end

function (f::CallWithMetadata)(args...)
    metadata(unwrap(f.f(map(unwrap, args)...)), metadata(f))
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
    ndim = :($length(($(indices...),)))

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

julia> @variables t z(t)[1:3] # also works for dependent variables
2-element Vector{Any}:
 t
  (z(t))[1:3]
```

A symbol or expression that represents an array can be turned into an array of
symbols or expressions using the `scalarize` function.

```julia
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

```julia
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
    esc(_parse_vars(:variables, Real, xs))
end

const _fail = Dict()

_getname(x, _) = nameof(x)
_getname(x::Symbol, _) = x
function _getname(x::Symbolic, val)
    issym(x) && return nameof(x)
    if iscall(x) && issym(operation(x))
        return nameof(operation(x))
    end
    if !hasmetadata(x, Symbolics.GetindexParent) && iscall(x) && operation(x) == getindex
        return _getname(arguments(x)[1], val)
    end
    ss = getsource(x, nothing)
    if ss === nothing
        ss = getsource(getparent(x), val)
    end
    ss === _fail && throw(ArgumentError("Variable $x doesn't have a source defined."))
    ss[2]
end

getsource(x, val=_fail) = getmetadata(unwrap(x), VariableSource, val)

SymbolicIndexingInterface.symbolic_type(::Type{<:Symbolics.Num}) = ScalarSymbolic()
SymbolicIndexingInterface.symbolic_type(::Type{<:Symbolics.Arr}) = ArraySymbolic()

SymbolicIndexingInterface.hasname(x::Union{Num,Arr}) = hasname(unwrap(x))

function SymbolicIndexingInterface.hasname(x::Symbolic)
    issym(x) || !iscall(x) || iscall(x) && (issym(operation(x)) || operation(x) == getindex)
end

SymbolicIndexingInterface.getname(x, val=_fail) = _getname(unwrap(x), val)

function SymbolicIndexingInterface.symbolic_evaluate(ex::Union{Num, Arr, Symbolic, Equation, Inequality}, d::Dict; kwargs...)
    fixpoint_sub(ex, d; kwargs...)
end

"""
    fixpoint_sub(expr, dict; operator = Nothing)

Given a symbolic expression, equation or inequality `expr` perform the substitutions in
`dict` recursively until the expression does not change. Substitutions that depend on one
another will thus be recursively expanded. For example,
`fixpoint_sub(x, Dict(x => y, y => 3))` will return `3`. The `operator` keyword can be
specified to prevent substitution of expressions inside operators of the given type.

See also: [`fast_substitute`](@ref).
"""
function fixpoint_sub(x, dict; operator = Nothing)
    y = fast_substitute(x, dict; operator)
    while !isequal(x, y)
        y = x
        x = fast_substitute(y, dict; operator)
    end

    return x
end

const Eq = Union{Equation, Inequality}
"""
    fast_substitute(expr, dict; operator = Nothing)

Given a symbolic expression, equation or inequality `expr` perform the substitutions in
`dict`. This only performs the substitutions once. For example,
`fast_substitute(x, Dict(x => y, y => 3))` will return `y`. The `operator` keyword can be
specified to prevent substitution of expressions inside operators of the given type.

See also: [`fixpoint_sub`](@ref).
"""
function fast_substitute(eq::Eq, subs; operator = Nothing)
    if eq isa Inequality
        Inequality(fast_substitute(eq.lhs, subs; operator),
            fast_substitute(eq.rhs, subs; operator),
            eq.relational_op)
    else
        Equation(fast_substitute(eq.lhs, subs; operator),
            fast_substitute(eq.rhs, subs; operator))
    end
end
function fast_substitute(eq::T, subs::Pair; operator = Nothing) where {T <: Eq}
    T(fast_substitute(eq.lhs, subs; operator), fast_substitute(eq.rhs, subs; operator))
end
function fast_substitute(eqs::AbstractArray, subs; operator = Nothing)
    fast_substitute.(eqs, (subs,); operator)
end
function fast_substitute(eqs::AbstractArray, subs::Pair; operator = Nothing)
    fast_substitute.(eqs, (subs,); operator)
end
for (exprType, subsType) in Iterators.product((Num, Symbolics.Arr), (Any, Pair))
    @eval function fast_substitute(expr::$exprType, subs::$subsType; operator = Nothing)
        fast_substitute(value(expr), subs; operator)
    end
end
function fast_substitute(expr, subs; operator = Nothing)
    if (_val = get(subs, expr, nothing)) !== nothing
        return _val
    end
    iscall(expr) || return expr
    op = fast_substitute(operation(expr), subs; operator)
    args = SymbolicUtils.arguments(expr)
    if !(op isa operator)
        canfold = Ref(!(op isa Symbolic))
        args = let canfold = canfold
            map(args) do x
                x′ = fast_substitute(x, subs; operator)
                canfold[] = canfold[] && !(x′ isa Symbolic)
                x′
            end
        end
        canfold[] && return op(args...)
    end
    maketerm(typeof(expr),
        op,
        args,
        symtype(expr),
        metadata(expr))
end
function fast_substitute(expr, pair::Pair; operator = Nothing)
    a, b = pair
    isequal(expr, a) && return b
    if a isa AbstractArray
        for (ai, bi) in zip(a, b)
            expr = fast_substitute(expr, ai => bi; operator)
        end
    end
    iscall(expr) || return expr
    op = fast_substitute(operation(expr), pair; operator)
    args = SymbolicUtils.arguments(expr)
    if !(op isa operator)
        canfold = Ref(!(op isa Symbolic))
        args = let canfold = canfold
            map(args) do x
                x′ = fast_substitute(x, pair; operator)
                canfold[] = canfold[] && !(x′ isa Symbolic)
                x′
            end
        end
        canfold[] && return op(args...)
    end
    maketerm(typeof(expr),
        op,
        args,
        symtype(expr),
        metadata(expr))
end

function getparent(x, val=_fail)
    maybe_parent = getmetadata(x, Symbolics.GetindexParent, nothing)
    if maybe_parent !== nothing
        return maybe_parent
    else
        if iscall(x) && operation(x) === getindex
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
    v = Sym{T}(name_ij)
    if T <: FnType
        v = CallWithMetadata(v)
    end
    Num(setmetadata(v, VariableSource, (:variables, name_ij)))
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

function rename_metadata(from, to, name)
    if hasmetadata(from, VariableSource)
        s = getmetadata(from, VariableSource)
        to = setmetadata(to, VariableSource, (s[1], name))
    end
    if hasmetadata(from, GetindexParent)
        s = getmetadata(from, GetindexParent)
        to = setmetadata(to, GetindexParent, rename(s, name))
    end
    return to
end

rename(x::Union{Num, Arr}, name) = wrap(rename(unwrap(x), name))
function rename(x::ArrayOp, name)
    t = x.term
    args = arguments(t)
    # Hack:
    @assert operation(t) === (map) && (args[1] isa CallWith || args[1] == CallWithMetadata)
    rn = rename(args[2], name)

    xx = metadata(operation(t)(args[1], rn), metadata(x))
    rename_getindex_source(rename_metadata(x, xx, name))
end

function rename(x::CallWithMetadata, name)
    rename_metadata(x, CallWithMetadata(rename(x.f, name), x.metadata), name)
end

function rename(x::Symbolic, name)
    if issym(x)
        xx = @set! x.name = name
        xx = rename_metadata(x, xx, name)
        symtype(xx) <: AbstractArray ? rename_getindex_source(xx) : xx
    elseif iscall(x) && operation(x) === getindex
        rename(arguments(x)[1], name)[arguments(x)[2:end]...]
    elseif iscall(x) && symtype(operation(x)) <: FnType || operation(x) isa CallWithMetadata
        xx = @set x.f = rename(operation(x), name)
        @set! xx.hash = Ref{UInt}(0)
        return rename_metadata(x, xx, name)
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
