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
    get_variables(e, varlist = nothing; sort::Bool = false)

Return a vector of variables appearing in `e`, optionally restricting to variables in `varlist`.

Note that the returned variables are not wrapped in the `Num` type.

# Examples
```jldoctest
julia> @variables t x y z(t);

julia> Symbolics.get_variables(x + y + sin(z); sort = true)
3-element Vector{SymbolicUtils.BasicSymbolic}:
 x
 y
 z(t)

julia> Symbolics.get_variables(x - y; sort = true)
2-element Vector{SymbolicUtils.BasicSymbolic}:
 x
 y
```
"""
function get_variables(e::Num, varlist = nothing; sort::Bool = false)
    get_variables(value(e), varlist; sort)
end
function get_variables(e, varlist = nothing; sort::Bool = false)
    vars = Vector{BasicSymbolic}()
    get_variables!(vars, e, varlist)
    if sort
        sort!(vars; by = SymbolicUtils.get_degrees)
    end
    vars
end

get_variables!(vars, e::Num, varlist=nothing) = get_variables!(vars, value(e), varlist)
get_variables!(vars, e, varlist=nothing) = vars

function is_singleton(e)
    if iscall(e)
        op = operation(e)
        op === getindex && return true
        iscall(op) && return is_singleton(op) # recurse to reach getindex for array element variables
        return issym(op) && !hasmetadata(e, CallWithParent)
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
        get_variables!(vars, operation(e), varlist)
        foreach(x -> get_variables!(vars, x, varlist), arguments(e))
    end
    return (vars isa AbstractVector) ? unique!(vars) : vars
end

function get_variables!(vars, e::Equation, varlist=nothing)
  get_variables!(vars, e.rhs, varlist)
end

# Sym / Term --> Symbol
Base.Symbol(x::Union{Num,Symbolic}) = tosymbol(x)
tosymbol(t::Num; kwargs...) = tosymbol(value(t); kwargs...)

"""
    diff2term(x, x_metadata::Dict{Datatype, Any}) -> Symbolic

Convert a differential variable to a `Term`. Note that it only takes a `Term`
not a `Num`.
Any upstream metadata can be passed via `x_metadata`

```jldoctest
julia> @variables x t u(x, t) z(t)[1:2]; Dt = Differential(t); Dx = Differential(x);

julia> Symbolics.diff2term(Symbolics.value(Dx(Dt(u))))
uˍtx(x, t)

julia> Symbolics.diff2term(Symbolics.value(Dt(z[1])))
(zˍt(t))[1]
```
"""
function diff2term(O, O_metadata::Union{Dict, Nothing, Base.ImmutableDict}=nothing)
    iscall(O) || return O

    inner = O
    op = identity
    while is_derivative(inner)
        op = op ∘ operation(inner)
        inner = arguments(inner)[1]
    end
    if iscall(inner) && operation(inner) == getindex
        return diff2term(op(arguments(inner)[1]), O_metadata)[arguments(inner)[2:end]...]
    end
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
        return maketerm(typeof(O), TermInterface.head(O), map(diff2term, children(O)),
                        O_metadata isa Nothing ?
            metadata(O) : Base.ImmutableDict(metadata(O)..., O_metadata...))
    else
        oldop = operation(O)
        opname = if issym(oldop)
            string(nameof(oldop))
        elseif iscall(oldop) && operation(oldop) === getindex
            string(nameof(arguments(oldop)[1]))
        elseif oldop isa Function
            return nothing
        else
            error("diff2term case not handled: $oldop")
        end
        newname = occursin(d_separator, opname) ? Symbol(opname, ds) : Symbol(opname, d_separator, ds)
        return setname(maketerm(typeof(O), rename(oldop, newname), children(O), O_metadata isa Nothing ?
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

```jldoctest
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
    elseif iscall(t)
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
    if issym(x) || x isa CallWithMetadata
        (x, i)
    elseif iscall(x)
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

```jldoctest
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

```jldoctest
julia> @variables a x y;

julia> Symbolics.coeff(2a, x)
0

julia> Symbolics.coeff(3x + 2y, y)
2

julia> Symbolics.coeff(x^2 + y, x^2)
1

julia> Symbolics.coeff(2*x*y + y, x*y)
2
```
"""
function coeff(p, sym=nothing)
    # if `sym` is a product, iteratively compute the coefficient w.r.t. each term in `sym`
    if iscall(value(sym)) && operation(value(sym)) === (*)
        for t in arguments(value(sym))
            @assert !(t isa Number) "`coeff(p, sym)` does not allow `sym` containing numerical factors"
            p = coeff(p, t)
        end
        return p
    end
            
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
        args = arguments(p)
        coeffs = map(a->coeff(a, sym), args)
        if all(_iszero, coeffs)
            return 0
        else
            @views prod(Iterators.flatten((coeffs[findall(!_iszero, coeffs)], args[findall(_iszero, coeffs)])))
        end
    elseif isdiv(p)
        numerator, denominator = arguments(p)
        if !occursin(sym, denominator)
            coeff(numerator, sym) / denominator
        else
            throw(DomainError(p, "coeff on fractions is not yet implemented."))
        end
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

"""
    symbolic_to_float(x::Union{Num, BasicSymbolic})::Union{AbstractFloat, BasicSymbolic}

If the symbolic value is exactly equal to a number, converts the symbolic value
to a floating point number. Otherwise retains the symbolic value.

## Examples

```julia
symbolic_to_float((1//2 * x)/x) # 0.5
symbolic_to_float((1/2 * x)/x) # 0.5
symbolic_to_float((1//2)*√(279//4)) # 4.175823272122517
```
"""
function symbolic_to_float end
symbolic_to_float(x::Num) = symbolic_to_float(unwrap(x))
symbolic_to_float(x::Number) = x
function symbolic_to_float(x::SymbolicUtils.BasicSymbolic)
    substitute(x,Dict())
end

"""
    numerator(x)

Return the numerator of the symbolic expression `x`.

Examples
========
```julia-repl
julia> numerator(x/y)
x
```
"""
function Base.numerator(x::Union{Num, Symbolic})
    x = unwrap(x)
    if iscall(x) && operation(x) == /
        x = arguments(x)[1] # get numerator
    end
    return wrap(x)
end

"""
    denominator(x)

Return the denominator of the symbolic expression `x`.

Examples
========
```julia-repl
julia> denominator(x/y)
y
```
"""
function Base.denominator(x::Union{Num, Symbolic})
    x = unwrap(x)
    if iscall(x) && operation(x) == /
        x = arguments(x)[2] # get denominator
    else
        x = 1
    end
    return wrap(x)
end

"""
    arguments(x, op::Function)

Get the arguments of the symbolic expression `x` with respect to the operation or function `op`.
"""
function arguments(x, op::Function)
    x = unwrap(x)
    if iscall(x) && operation(x) == op
        args = [arguments(arg, op) for arg in arguments(x)] # recurse into each argument and obtain its factors
        args = reduce(vcat, args) # concatenate array of arrays into one array
    else
        args = [wrap(x)] # base case
    end
    return args
end

"""
    terms(x)

Get the terms of the symbolic expression `x`.

Examples
========
```julia-repl
julia> terms(-x + y - z)
3-element Vector{Num}:
 -z
  y
 -x
```
"""
terms(x) = arguments(x, +)

"""
    factors(x)

Get the factors of the symbolic expression `x`.

Examples
========
```julia-repl
julia> factors(2 * x * y)
3-element Vector{Num}:
 2
 y
 x
```
"""
factors(x) = arguments(x, *)

"""
    evaluate(eq::Equation, subs)
    evaluate(ineq::Inequality, subs)

Evaluate the equation `eq` or inequality `ineq`. `subs` is a dictionary of variable to numerical value substitutions. 
If both sides of the equation or inequality are numeric, then the result is a boolean. 

# Examples
```julia-repl
julia> @variables x y
julia> eq = x ~ y
julia> evaluate(eq, Dict(x => 1, y => 1))
true

julia> ltr = x ≲ y
julia> evaluate(ltr, Dict(x => 1, y => 2))
true

julia> gtr = x ≳ y
julia> evaluate(gtr, Dict(x => 1, y => 2))
false
```
"""
function evaluate end

function evaluate(eq::Equation, subs)
    lhs = fast_substitute(eq.lhs, subs)
    rhs = fast_substitute(eq.rhs, subs)
    return isequal(lhs, rhs)
end

function evaluate(ineq::Inequality, subs)
    lhs = fast_substitute(ineq.lhs, subs)
    rhs = fast_substitute(ineq.rhs, subs)
    if (ineq.relational_op == geq)
        return isless(rhs, lhs)
    elseif (ineq.relational_op == leq)
        return isless(lhs, rhs)
    else
        throw(ArgumentError("Inequality $ineq not supported"))
    end
end 


