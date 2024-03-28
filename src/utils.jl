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
    get_variables(O) -> Vector{BasicSymbolic}

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

function is_singleton(e)
    if istree(e)
        op = operation(e)
        op === getindex && return true
        istree(op) && return is_singleton(op) # recurse to reach getindex for array element variables
        return issym(op)
    else
        return issym(e)
    end
end

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
tosymbol(t::Num; kwargs...) = tosymbol(value(t); kwargs...)

"""
    diff2term(x, x_metadata::Dict{Datatype, Any}) -> Symbolic

Convert a differential variable to a `Term`. Note that it only takes a `Term`
not a `Num`.
Any upstream metadata can be passed via `x_metadata`

```julia
julia> @variables x t u(x, t) z(t)[1:2]; Dt = Differential(t); Dx = Differential(x);

julia> Symbolics.diff2term(Symbolics.value(Dx(Dt(u))))
uˍtx(x, t)

julia> Symbolics.diff2term(Symbolics.value(Dt(z[1])))
var"z(t)[1]ˍt"
```
"""
function diff2term(O, O_metadata::Union{Dict, Nothing, Base.ImmutableDict}=nothing)
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
        return similarterm(O, operation(O), map(diff2term, arguments(O)), metadata = O_metadata isa Nothing ?
            metadata(O) : Base.ImmutableDict(metadata(O)..., O_metadata...))
    else
        oldop = operation(O)
        if issym(oldop)
            opname = string(nameof(oldop))
        elseif istree(oldop) && operation(oldop) === getindex
            opname = string(nameof(arguments(oldop)[1]))
            args = arguments(O)
        elseif oldop == getindex
            args = arguments(O)
            opname = string(tosymbol(args[1]), "[", map(tosymbol, args[2:end])..., "]")
            return Sym{symtype(O)}(Symbol(opname, d_separator, ds))
        end
        newname = occursin(d_separator, opname) ? Symbol(opname, ds) : Symbol(opname, d_separator, ds)
        return setname(similarterm(O, rename(oldop, newname), arguments(O), metadata = O_metadata isa Nothing ?
            metadata(O) : Base.ImmutableDict(metadata(O)..., O_metadata...)), newname)
    end
end

setname(v, name) = setmetadata(v, Symbolics.VariableSource, (:variables, name))

"""
    tosymbol(x::Union{Num,Symbolic}; states=nothing, escape=true) -> Symbol

Convert `x` to a symbol. `states` are the states of a system, and `escape`
means if the target has escapes like `val"y(t)"`. If `escape` is false, then
it will only output `y` instead of `y(t)`.

# Examples

```julia
julia> @variables t z(t)
2-element Vector{Num}:
    t
 z(t)

julia> Symbolics.tosymbol(z)
Symbol("z(t)")

julia> Symbolics.tosymbol(z; escape=false)
:z
```
"""
function tosymbol(t; states=nothing, escape=true)
    if issym(t)
        return nameof(t)
    elseif istree(t)
        if issym(operation(t))
            if states !== nothing && !(t in states)
                return nameof(operation(t))
            end
            op = nameof(operation(t))
            args = arguments(t)
        elseif operation(t) isa Differential
            term = diff2term(t)
            if issym(term)
                return nameof(term)
            else
                op = Symbol(operation(term))
                args = arguments(term)
            end
        else
            op = Symbol(repr(operation(t)))
            args = arguments(t)
        end

        return escape ? Symbol(op, "(", join(args, ", "), ")") : op
    else
        return t
    end
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

"""
    degree(p, sym=nothing)

Extract the degree of `p` with respect to `sym`.

# Examples

```julia
julia> @variables x;

julia> Symbolics.degree(x^0)
0

julia> Symbolics.degree(x)
1

julia> Symbolics.degree(x^2)
2
```
"""
function degree(p, sym=nothing)
    p = value(p)
    sym = value(sym)
    if p isa Number
        return 0
    end
    if isequal(p, sym)
        return 1
    end
    if isterm(p)
        if sym === nothing
            return 1
        else
            return Int(isequal(p, sym))
        end
    elseif ismul(p)
        return sum(degree(k^v, sym) for (k, v) in zip(keys(p.dict), values(p.dict)))
    elseif isadd(p)
        return maximum(degree(key, sym) for key in keys(p.dict))
    elseif ispow(p)
        return p.exp * degree(p.base, sym)
    elseif isdiv(p)
        return degree(p.num, sym) - degree(p.den, sym)
    elseif issym(p)
        if sym === nothing
            return 1
        else
            return Int(isequal(p, sym))
        end
    end
    throw(DomainError(p, "Datatype $(typeof(p)) not accepted."))
end

"""
    coeff(p, sym=nothing)

Extract the coefficient of `p` with respect to `sym`.
Note that `p` might need to be expanded and/or simplified with `expand` and/or `simplify`.

# Examples

```julia
julia> @variables a x y;

julia> Symbolics.coeff(2a, x)
0

julia> Symbolics.coeff(3x + 2y, y)
2

julia> Symbolics.coeff(x^2 + y, x^2)
1
```
"""
function coeff(p, sym=nothing)
    p, sym = value(p), value(sym)

    if isequal(sym, 1)
        sym = nothing
    end

    if issym(p) || isterm(p)
        sym === nothing ? 0 : Int(isequal(p, sym))
    elseif ispow(p)
        sym === nothing ? 0 : Int(isequal(p, sym))
    elseif isadd(p)
        if sym===nothing
            p.coeff
        else
            sum(coeff(k, sym) * v for (k, v) in p.dict)
        end
    elseif ismul(p)
        args = unsorted_arguments(p)
        coeffs = map(a->coeff(a, sym), args)
        if all(iszero, coeffs)
            return 0
        else
            @views prod(Iterators.flatten((coeffs[findall(!iszero, coeffs)], args[findall(iszero, coeffs)])))
        end
    elseif isdiv(p)
        throw(DomainError(p, "coeff on expressions with division is not yet implemented."))
    else
        p isa Number && return sym === nothing ? p : 0
        p isa Symbolic && return coeff(p, sym)
        throw(DomainError(p, "Datatype $(typeof(p)) not accepted."))
    end
end

### Nums <--> Polys

const DP = DynamicPolynomials
# extracting underlying polynomial and coefficient type from Polyforms
underlyingpoly(x::Number) = x
underlyingpoly(pf::PolyForm) = pf.p
coefftype(x::Number) = typeof(x)
coefftype(pf::PolyForm) = DP.coefficient_type(underlyingpoly(pf))

#=
Converts an array of symbolic polynomials
into an array of DynamicPolynomials.Polynomials
=#
function symbol_to_poly(sympolys::AbstractArray)
    @assert !isempty(sympolys) "Empty input."

    # standardize input
    stdsympolys = map(unwrap, sympolys)
    sort!(stdsympolys, lt=(<ₑ))

    pvar2sym = Bijections.Bijection{Any,Any}()
    sym2term = Dict{BasicSymbolic,Any}()
    polyforms = map(f -> PolyForm(f, pvar2sym, sym2term), stdsympolys)

    # Discover common coefficient type
    commontype = mapreduce(coefftype, promote_type, polyforms, init=Int)
    @assert commontype <: Union{Integer,Rational} "Only integer and rational coefficients are supported as input."

    # Convert all to DP.Polynomial, so that coefficients are of same type,
    # and constants are treated as polynomials
    # We also need this because Groebner does not support abstract types as input
    polynoms = Vector{DP.Polynomial{DP.Commutative{DP.CreationOrder},DP.Graded{DP.LexOrder},commontype}}(undef, length(sympolys))
    for (i, pf) in enumerate(polyforms)
        polynoms[i] = underlyingpoly(pf)
    end

    polynoms, pvar2sym, sym2term
end

#=
Converts an array of AbstractPolynomialLike`s into an array of
symbolic expressions mapping variables w.r.t pvar2sym
=#
function poly_to_symbol(polys, pvar2sym, sym2term, ::Type{T}) where {T}
    map(f -> PolyForm{T}(f, pvar2sym, sym2term), polys)
end
