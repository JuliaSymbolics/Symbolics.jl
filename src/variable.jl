using SymbolicUtils: FnType, Sym
using Setfield

const IndexMap = Dict{Char,Char}(
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

"""
$(TYPEDEF)

A named variable which represents a numerical value. The variable is uniquely
identified by its `name`, and all variables with the same `name` are treated
as equal.

# Fields
$(FIELDS)

For example, the following code defines an independent variable `t`, a parameter
`α`, a function parameter `σ`, a variable `x`, which depends on `t`, a variable
`y` with no dependents, a variable `z`, which depends on `t`, `α`, and `x(t)`
and parameters `β₁` and `β₂`.


```julia
σ = Num(Variable{Symbolics.FnType{Tuple{Any},Real}}(:σ)) # left uncalled, since it is used as a function
w = Num(Variable{Symbolics.FnType{Tuple{Any},Real}}(:w)) # unknown, left uncalled
x = Num(Variable{Symbolics.FnType{Tuple{Any},Real}}(:x))(t)  # unknown, depends on `t`
y = Num(Variable(:y))   # unknown, no dependents
z = Num(Variable{Symbolics.FnType{NTuple{3,Any},Real}}(:z))(t, α, x)  # unknown, multiple arguments
β₁ = Num(Variable(:β, 1)) # with index 1
β₂ = Num(Variable(:β, 2)) # with index 2

expr = β₁ * x + y^α + σ(3) * (z - t) - β₂ * w(t - 1)
```
"""
struct Variable{T} <: Function # backward compat
    """The variable's unique name."""
    name::Symbol
    Variable(name) = Sym{Real}(name)
    Variable{T}(name) where T = Sym{T}(name)
    function Variable{T}(name, indices...) where T
        var_name = Symbol("$(name)$(join(map_subscripts.(indices), "ˏ"))")
        Sym{T}(var_name)
    end
end

function Variable(name, indices...)
    var_name = Symbol("$(name)$(join(map_subscripts.(indices), "ˏ"))")
    Variable(var_name)
end

# TODO: move this to Symutils
function Sym{T}(name, i, indices...) where T
    var_name = Symbol("$(name)$(join(map_subscripts.((i, indices...,)), "ˏ"))")
    Sym{T}(var_name)
end

function map_subscripts(indices)
    str = string(indices)
    join(IndexMap[c] for c in str)
end

rename(x::Sym,name) = @set! x.name = name
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

        # x [connect = flow; unit = u"m^3/s"]
        if isoption(nv)
            options = nv.args
            cursor += 1
        end

        iscall = Meta.isexpr(v, :call)
        isarray = Meta.isexpr(v, :ref)
        issym  = v isa Symbol
        @assert iscall || isarray || issym "@$macroname expects a tuple of expressions or an expression of a tuple (`@$macroname x y z(t) v[1:3] w[1:2,1:4]` or `@$macroname x y z(t) v[1:3] w[1:2,1:4] k=1.0`)"

        if iscall
            var_name, expr = construct_vars(v.args[1], type, v.args[2:end], val, options, transform)
        else
            var_name, expr = construct_vars(v, type, nothing, val, options, transform)
        end

        push!(var_names, var_name)
        push!(ex.args, expr)
    end
    rhs = build_expr(:tuple, var_names)
    push!(ex.args, :(($(var_names...),) = $rhs))
    return ex
end

function construct_vars(v, type, call_args, val, prop, transform)
    issym  = v isa Symbol
    isarray = isa(v, Expr) && v.head == :ref
    if isarray
        var_name = v.args[1]
        indices = v.args[2:end]
        expr = _construct_array_vars(var_name, type, call_args, val, prop, indices...)
    else
        var_name = v
        expr = construct_var(var_name, type, call_args, val, prop)
    end
    var_name, :($var_name = $transform($expr))
end

function option_to_metadata_type(::Val{opt}) where {opt}
    throw(Base.Meta.ParseError("unknown property type $opt"))
end

function setprops_expr(expr, props)
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


function construct_var(var_name, type, call_args, val, prop)
    expr = if call_args === nothing
        :($Num($Sym{$type}($(Meta.quot(var_name)))))
    elseif !isempty(call_args) && call_args[end] == :..
        :($Num($Sym{$FnType{Tuple, $type}}($(Meta.quot(var_name))))) # XXX: using Num as output
    else
        :($Num($Sym{$FnType{NTuple{$(length(call_args)), Any}, $type}}($(Meta.quot(var_name)))($(map(x->:($value($x)), call_args)...))))
    end

    if val !== nothing
        expr = :($setmetadata($expr, $VariableDefaultValue, $val))
    end

    return setprops_expr(expr, prop)
end

function construct_var(var_name, type, call_args, val, prop, ind)
    # TODO: just use Sym here
    expr = if call_args === nothing
        :($Num($Sym{$type}($(Meta.quot(var_name)), $ind...)))
    elseif !isempty(call_args) && call_args[end] == :..
        :($Num($Sym{$FnType{Tuple{Any}, $type}}($(Meta.quot(var_name)), $ind...))) # XXX: using Num as output
    else
        :($Num($Sym{$FnType{NTuple{$(length(call_args)), Any}, $type}}($(Meta.quot(var_name)), $ind...)($(map(x->:($value($x)), call_args)...))))
    end
    if val !== nothing
        expr = :($setmetadata($expr, $VariableDefaultValue, $val isa AbstractArray ? $val[$ind...] : $val))
    end

    return setprops_expr(expr, prop)
end

function _construct_array_vars(var_name, type, call_args, val, prop, indices...)
    :(map(Iterators.product($(indices...))) do ind
        $(construct_var(var_name, type, call_args, val, prop, :ind))
    end)
end


"""

Define one or more unknown variables.

```julia
@variables t α σ(..) β[1:2]
@variables w(..) x(t) y z(t, α, x)

expr = β₁* x + y^α + σ(3) * (z - t) - β₂ * w(t - 1)
```

`(..)` signifies that the value should be left uncalled.

Sometimes it is convenient to define arrays of variables to model things like `x₁,…,x₃`.
The `@variables` macro supports this with the following syntax:

```julia
@variables x[1:3];
x

3-element Vector{Num}:
 x₁
 x₂
 x₃

# support for arbitrary ranges and tensors
@variables y[2:3,1:5:6];
y

2×2 Matrix{Num}:
 y₂ˏ₁  y₂ˏ₆
 y₃ˏ₁  y₃ˏ₆

# also works for dependent variables
@variables t z[1:3](t);
z

3-element Array{Num,1}:
 z₁(t)
 z₂(t)
 z₃(t)
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
