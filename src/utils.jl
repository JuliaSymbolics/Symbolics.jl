isblock(x) = length(x) == 1 && x[1] isa Expr && x[1].head == :block
function flatten_expr!(x)
    isblock(x) || return x
    x = MacroTools.striplines(x[1])
    filter!(z -> z isa Symbol || z.head != :line, x.args)
    xs = []
    for ex in x.args
        if Meta.isexpr(ex, :tuple)
            append!(xs, ex.args)
        else
            push!(xs, ex)
        end
    end
    xs
end
function build_expr(head::Symbol, args)
    ex = Expr(head)
    append!(ex.args, args)
    ex
end

"""
    get_variables(O) -> Vector{Union{Sym, Term}}

Returns the variables in the expression. Note that the returned variables are
not wrapped in the `Num` type.

# Examples
```julia
julia> @variables t x y z(t)
4-element Vector{Num}:
    t
    x
    y
 z(t)

julia> ex = x + y + sin(z)
(x + y) + sin(z(t))

julia> Symbolics.get_variables(ex)
3-element Vector{Any}:
 x
 y
 z(t)
```
"""
get_variables(e::Num, varlist=nothing) = get_variables(value(e), varlist)
get_variables!(vars, e::Num, varlist=nothing) = get_variables!(vars, value(e), varlist)
get_variables!(vars, e, varlist=nothing) = vars

function is_singleton(e::Term)
    op = operation(e)
    op === getindex && return true
    op isa Term && return is_singleton(op) # recurse to reach getindex for array element variables
    op isa Sym
end

is_singleton(e::Sym) = true
is_singleton(e) = false

get_variables!(vars, e::Number, varlist=nothing) = vars

function get_variables!(vars, e::Symbolic, varlist=nothing)
    if is_singleton(e)
        if isnothing(varlist) || any(isequal(e), varlist)
            push!(vars, e)
        end
    else
        foreach(x -> get_variables!(vars, x, varlist), arguments(e))
    end
    return (vars isa AbstractVector) ? unique!(vars) : vars
end

function get_variables!(vars, e::Equation, varlist=nothing)
  get_variables!(vars, e.rhs, varlist)
end

get_variables(e, varlist=nothing) = get_variables!([], e, varlist)

# Sym / Term --> Symbol
Base.Symbol(x::Union{Num,Symbolic}) = tosymbol(x)
tosymbol(x; kwargs...) = x
tosymbol(x::Sym; kwargs...) = nameof(x)
tosymbol(t::Num; kwargs...) = tosymbol(value(t); kwargs...)

"""
    diff2term(x::Term) -> Symbolic
    diff2term(x) -> x

Convert a differential variable to a `Term`. Note that it only takes a `Term`
not a `Num`.

```julia
julia> @variables x t u(x, t) z(t)[1:2]; Dt = Differential(t); Dx = Differential(x);

julia> Symbolics.diff2term(Symbolics.value(Dx(Dt(u))))
uˍtx(x, t)

julia> Symbolics.diff2term(Symbolics.value(Dt(z[1])))
var"z(t)[1]ˍt"
```
"""
function diff2term(O)
    istree(O) || return O
    if is_derivative(O)
        ds = ""
        while is_derivative(O)
            ds = string(operation(O).x) * ds
            O = arguments(O)[1]
        end
    else
        ds = nothing
    end
    d_separator = 'ˍ'

    if ds === nothing
        return similarterm(O, operation(O), map(diff2term, arguments(O)), metadata=metadata(O))
    else
        oldop = operation(O)
        if oldop isa Sym
            opname = string(nameof(oldop))
            args = arguments(O)
        elseif oldop isa Term && operation(oldop) === getindex
            opname = string(nameof(arguments(oldop)[1]))
            args = arguments(O)
        elseif oldop == getindex
            args = arguments(O)
            opname = string(tosymbol(args[1]), "[", map(tosymbol, args[2:end])..., "]")
            return Sym{symtype(O)}(Symbol(opname, d_separator, ds))
        end
        newname = occursin(d_separator, opname) ? Symbol(opname, ds) : Symbol(opname, d_separator, ds)
        return setname(similarterm(O, rename(oldop, newname), arguments(O), metadata=metadata(O)), newname)
    end
end

setname(v, name) = setmetadata(v, Symbolics.VariableSource, (:variables, name))

"""
    tosymbol(x::Union{Num,Symbolic}; states=nothing, escape=true) -> Symbol

Convert `x` to a symbol. `states` are the states of a system, and `escape`
means if the target has escapes like `val"y(t)"`. If `escape` is false then
it will only output `y` instead of `y(t)`.

# Examples

```julia
julia> @variables t z(t)
2-element Vector{Num}:
    t
 z(t)

julia> Symbolics.tosymbol(z)
Symbol("z(t)")

julia>  Symbolics.tosymbol(z; escape=false)
:z
```
"""
function tosymbol(t::Term; states=nothing, escape=true)
    if operation(t) isa Sym
        if states !== nothing && !(t in states)
            return nameof(operation(t))
        end
        op = nameof(operation(t))
        args = arguments(t)
    elseif operation(t) isa Differential
        term = diff2term(t)
        if issym(term)
            return nameof(term)
        end
        op = Symbol(operation(term))
        args = arguments(term)
    else
        op = Symbol(repr(operation(t)))
        args = arguments(t)
    end

    return escape ? Symbol(op, "(", join(args, ", "), ")") : op
end

function lower_varname(var::Symbolic, idv, order)
    order == 0 && return var
    D = Differential(idv)
    for _ in 1:order
        var = D(var)
    end
    return diff2term(var)
end

### OOPS

struct Unknown end

macro oops(ex)
    quote
        tmp = $(esc(ex))
        if tmp === Unknown()
            return Unknown()
        else
            tmp
        end
    end
end

function makesubscripts(n)
    set = 'i':'z'
    m = length(set)
    map(1:n) do i
        repeats = ceil(Int, i / m)
        c = set[(i-1) % m + 1]
        Sym{Int}(Symbol(join([c for _ in 1:repeats], "")))
    end
end

function var_from_nested_derivative(x,i=0)
    x = unwrap(x)
    if issym(x)
        (x, i)
    elseif istree(x)
        operation(x) isa Differential ?
            var_from_nested_derivative(first(arguments(x)), i + 1) : (x, i)
    else
        error("Not a well formed derivative expression $x")
    end
end

degree(p::Union{Term,Sym}, sym=nothing) = sym === nothing ? 1 : Int(isequal(p, sym))
degree(p::Add, sym=nothing) = maximum(degree.(keys(p.dict), sym))
degree(p::Mul, sym=nothing) = sum(degree(k, sym) * v for (k, v) in p.dict)
degree(p::Pow, sym=nothing) = p.exp * degree(p.base, sym)

"""
    degree(p, sym=nothing)

Extract the degree of `p` with respect to `sym`.
"""
function degree(p, sym=nothing)
    p, sym = value(p), value(sym)
    p isa Number && return 0
    isequal(p, sym) && return 1
    p isa Symbolic && return degree(p, sym)
    throw(DomainError(p, "Datatype $(typeof(p)) not accepted."))
end

coeff(p::Union{Term,Sym}, sym=nothing) = sym === nothing ? 0 : Int(isequal(p, sym))
coeff(p::Pow, sym=nothing) = sym === nothing ? 0 : Int(isequal(p, sym))
coeff(p::Add, sym=nothing) = sum(coeff.(arguments(p), sym))
function coeff(p::Mul, sym=nothing)
    args = arguments(p)
    I = findall(a -> !isequal(a, sym), args)
    length(I) == length(args) ? 0 : prod(args[I])
end

"""
    coeff(p, sym=nothing)

Extract the coefficient of `p` with respect to `sym`.
Note that `p` might need to be expanded and/or simplified with `expand` and/or `simplify`.
"""
function coeff(p, sym=nothing)
    p, sym = value(p), value(sym)
    p isa Number && return sym === nothing ? p : 0
    p isa Symbolic && return coeff(p, sym)
    throw(DomainError(p, "Datatype $(typeof(p)) not accepted."))
end
