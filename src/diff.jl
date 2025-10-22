"""
$(TYPEDEF)

Represents a differential operator.

# Fields
$(FIELDS)

# Examples

```jldoctest
julia> using Symbolics

julia> @variables x y;

julia> D = Differential(x)
(D'~x)

julia> D(y) # Differentiate y wrt. x
(D'~x)(y)

julia> Dx = Differential(x) * Differential(y) # d^2/dxy operator
(D'~x(t)) ∘ (D'~y(t))

julia> D3 = Differential(x)^3 # 3rd order differential operator
(D'~x(t)) ∘ (D'~x(t)) ∘ (D'~x(t))
```
"""
struct Differential <: Operator
    """The variable or expression to differentiate with respect to."""
    x::BasicSymbolic{VartypeT}
    """The derivative order. Can be rational for fractional derivatives."""
    order::Union{Int, Rational{Int}}
    function Differential(x::BasicSymbolic{VartypeT}, order = 1)
        @assert order > 0 "Derivative order must be positive"
        @match x begin
            BSImpl.Const(;) => throw(ArgumentError("Cannot take derivative with respect to constant."))
            _ => new(x, order)
        end
    end
    Differential(x::Union{Num, Arr}, order = 1) = Differential(unwrap(x), order)
    Differential(::CallAndWrap, order = 1) = throw(ArgumentError("Cannot take derivative with respect to a symbolic function."))
    Differential(::Union{AbstractFloat, Integer}) = error("D(::Number) is not a valid derivative. Derivatives must be taken w.r.t. symbolic variables.")
end
function (D::Differential)(x::BasicSymbolic{VartypeT})
    @match x begin
        BSImpl.Term(; f, args) && if f isa Differential && isequal(f.x, D.x) end => begin
            return Differential(D.x, D.order + f.order)(args[1])
        end
        _ => return BSImpl.Term{VartypeT}(D, SArgsT((x,)); type = symtype(x), shape = shape(x))
    end
end

(D::Differential)(::Union{AbstractFloat, Integer}) = Num(COMMON_ZERO)
(D::Differential)(x::Num) = Num(D(unwrap(x)))
(D::Differential)(x::Arr{T, N}) where {T, N} = Arr{T, N}(D(unwrap(x)))
(D::Differential)(x::Complex{Num}) = Complex{Num}(Num(D(unwrap(real(x)))), Num(D(unwrap(imag(x)))))
SymbolicUtils.isbinop(f::Differential) = false

function (s::SymbolicUtils.Substituter)(x::Differential)
    Differential(s(x.x), x.order)
end

function SymbolicUtils.operator_to_term(::Differential, ex::BasicSymbolic{VartypeT})
    return diff2term(ex)
end

function is_derivative(x::SymbolicT)
    @match x begin
        BSImpl.Term(; f) && if f isa Differential end => true
        _ => false
    end
end
is_derivative(_) = false

Base.:*(D1::ComposedFunction, D2::Differential) = D1 ∘ D2
Base.:*(D1::Differential, D2) = D1 ∘ D2
Base.:*(D1::Differential, D2::Differential) = D1 ∘ D2
function Base.:^(D::Differential, n::Integer)
    iszero(n) && return identity
    return Differential(D.x, D.order * n)
end

Base.show(io::IO, D::Differential) = print(io, "Differential(", D.x, ", ", D.order, ")")
Base.nameof(D::Differential) = :Differential

Base.:(==)(D1::Differential, D2::Differential) = isequal(D1.x, D2.x) && isequal(D1.order, D2.order)
Base.hash(D::Differential, u::UInt) = hash(D.order, hash(D.x, xor(u, 0xdddddddddddddddd)))

"""
    $(TYPEDSIGNATURES)

Clear caches of all cached functions involved in computing derivatives.
"""
function clear_derivative_caches!() # public
    SymbolicUtils.clear_cache!(occursin_info)
    SymbolicUtils.clear_cache!(recursive_hasoperator)
end

"""
    $(TYPEDSIGNATURES)

Toggle caching in derivative related functions.
"""
function toggle_derivative_caching!(value::Bool) # public
    SymbolicUtils.toggle_caching!(occursin_info, value)
    SymbolicUtils.toggle_caching!(recursive_hasoperator, value)
end

"""
    $(TYPEDSIGNATURES)

Return a collection of all functions involved in derivative computation
that are cached.
"""
function cached_derivative_functions() # public
    (occursin_info, recursive_hasoperator)
end

SymbolicUtils.@cache limit = 500_000 function occursin_info(x::BasicSymbolic{VartypeT}, expr::BasicSymbolic{VartypeT}, fail::Bool = true)::Bool
    _occursin_info(x, expr, fail)
end

occursin_info(x::BasicSymbolic{VartypeT}, expr, fail::Bool = true) = false

# Allow scalarized expressions
@inline function is_scalar_indexed(ex::BasicSymbolic{VartypeT})
    # any `BSImpl.Term(; f = getindex)` is scalar - slicing is an arrayop
    @match ex begin
        BSImpl.Term(; f) && if f === getindex end => true
        _ => false
    end
end

function _occursin_info(x::BasicSymbolic{VartypeT}, expr::BasicSymbolic{VartypeT}, fail::Bool = true)
    shexpr = shape(expr)
    if SymbolicUtils.is_array_shape(shexpr)
        fail && error("Differentiation with array expressions is not yet supported")
        return SymbolicUtils.query(isequal(x), expr)
    end

    iscall(expr) || return isequal(x, expr)
    isequal(x, expr) && return true

    isix = is_scalar_indexed(x)
    isie = is_scalar_indexed(expr)

    if isie
        isix && return false
        return SymbolicUtils.query(isequal(x), expr)
    end

    op = operation(expr)

    if op isa Integral
        # check if x occurs in limits
        domain = op.domain
        lower, upper = unwrap.(DomainSets.endpoints(domain.domain))
        (occursin_info(x, lower) || occursin_info(x, upper)) && return true

        # check if x is shadowed by integration variable in integrand
        isequal(domain.variables, x) && return false
    end

    predicate = let cond = op !== getindex, x = x
        function __predicate(a)
            occursin_info(x, a, cond)
        end
    end
    any(predicate, arguments(expr))
end

"""
    hasderiv(O)

Returns true if the expression or equation `O` contains [`Differential`](@ref) terms.
"""
hasderiv(O) = recursive_hasoperator(Differential, O)


_recursive_hasoperator(::Type{op}, eq::Equation) where {op} = recursive_hasoperator(op, eq.lhs) || recursive_hasoperator(op, eq.rhs)
_recursive_hasoperator(::Type{op}) where {op} = Base.Fix1(recursive_hasoperator, op) # curry version
_recursive_hasoperator(::Type{T}, ::T) where T = true


"""
    recursive_hasoperator(op, O)

An internal function that contains the logic for [`hasderiv`](@ref).
Return true if `O` contains a term with `Operator` `op`.
"""
SymbolicUtils.@cache function recursive_hasoperator(op::Any, O::Any)::Bool
    _recursive_hasoperator(op, O)
end

function _recursive_hasoperator(::Type{op}, O::SymbolicT) where {op}
    @match O begin
        BSImpl.Const(;) => false
        BSImpl.Sym(;) => false
        BSImpl.Term(; f, args) => f isa op || any(_recursive_hasoperator(op), args)
        BSImpl.AddMul(; dict) => any(_recursive_hasoperator(op), keys(dict))
        BSImpl.Div(; num, den) => recursive_hasoperator(op, num) || recursive_hasoperator(op, den)
        BSImpl.ArrayOp(; expr, term) => begin
            if term isa SymbolicT
                recursive_hasoperator(op, term) && return true
            end
            recursive_hasoperator(op, expr)
        end
    end
end
_recursive_hasoperator(::Type{op}, O) where {op} = false

struct DerivativeNotDefinedError <: Exception
    expr
    i::Int
end

function Base.showerror(io::IO, err::DerivativeNotDefinedError)
    op = operation(err.expr)
    nargs = length(arguments(err.expr))
    # `Markdown.parse` instead of `@md_str` to allow interpolating inside `literal` blocks
    # and code fences
    err_str = Markdown.parse("""
        Derivative of `$(err.expr)` with respect to its $(err.i)-th argument is not defined.
        Define a derivative by adding a method to `Symbolics.derivative`:

        ```julia
        function Symbolics.derivative(::typeof($op), args::NTuple{$nargs, Any}, ::Val{$(err.i)})
            # ...
        end
        ```

        Refer to the documentation for `Symbolics.derivative` and the
        "[Adding Analytical Derivatives](@ref)" section of the docs for further information.
        """)
    show(io, MIME"text/plain"(), err_str)
end

function symdiff_substitute_filter(ex::BasicSymbolic{T}) where {T}
    SymbolicUtils.default_substitute_filter(ex) || @match ex begin
        BSImpl.Term(; f) && if f isa Differential end => true
        _ => false
    end
end

"""
    $(TYPEDSIGNATURES)

Identical to `substitute` except it also substitutes inside `Differential` operator
applications.
"""
function substitute_in_deriv(ex, rules; kw...)
    substitute(ex, rules; kw..., filterer = symdiff_substitute_filter)
end

function chain_diff(D::Differential, arg::BasicSymbolic{VartypeT}, inner_args::SymbolicUtils.ROArgsT{VartypeT}; kw...)
    any(isequal(D.x), inner_args) && return D(arg)

    summed_args = SymbolicUtils.ArgsT{VartypeT}()
    sizehint!(summed_args, length(inner_args))
    for a in inner_args
        t1 = executediff(Differential(a), arg; kw...)
        t2 = executediff(D, a; kw...)
        push!(summed_args, t1 * t2)
    end
    return SymbolicUtils.add_worker(VartypeT, summed_args)
end

"""
    executediff(D, arg; simplify=false, occurrences=nothing)

Apply the passed Differential D on the passed argument.

This function differs to `expand_derivatives` in that in only expands the
passed differential and not any other Differentials it encounters.

# Arguments
- `D::Differential`: The differential to apply
- `arg::BasicSymbolic`: The symbolic expression to apply the differential on.
- `simplify::Bool=false`: Whether to simplify the resulting expression using
    [`SymbolicUtils.simplify`](@ref).
- `occurrences=nothing`: Information about the occurrences of the independent
    variable in the argument of the derivative. This is used internally for
    optimization purposes.
- `throw_no_derivative=false`: Whether to throw if a function with unknown
    derivative is encountered.
"""
function executediff(D::Differential, arg::BasicSymbolic{VartypeT}; simplify=false, throw_no_derivative=false)
    isinteger(D.order) || throw(ArgumentError("`executediff` requires integer derivative order."))
    order = floor(Int, D.order)
    if order > 1
        for _ in 1:order
            arg = executediff(Differential(D.x), arg; simplify, throw_no_derivative)
        end
        return arg
    end
    isequal(arg, D.x) && return COMMON_ONE
    occursin_info(D.x, arg) || return COMMON_ZERO

    # We can safely assume `arg` is scalar, else `occursin_info` would have errored.
    @match arg begin
        # Const case will never be reached because of `occursin_info`
        # if the sym were equal to `D.x` we wouldn't be here
        BSImpl.Sym(;) => return COMMON_ZERO
        BSImpl.Term(; f, args) => begin
            if f isa BasicSymbolic{VartypeT}
                # the only case where `f` is a symbolic is if this is a called symbolic
                # function or a dependent variable. In either case, we know it contains
                # `D.x` because of `occursin_info` and will just return `D(arg)`
                inner_args = arguments(arg)
                return chain_diff(D, arg, inner_args; simplify, throw_no_derivative)
            elseif f === getindex
                arr = arguments(arg)[1]
                inner_args = arguments(arguments(arg)[1])
                idx = SymbolicUtils.StableIndex(@views arguments(arg)[2:end])
                summed_args = SymbolicUtils.ArgsT{VartypeT}()
                sizehint!(summed_args, length(inner_args))
                # We know `D.x` is in `arg`, so the derivative is not identically zero.
                # `arg` cannot be `D.x` since, that would have also early exited. 
                for (i, a) in enumerate(inner_args)
                    der = derivative_idx(arr, i)
                    if isequal(a, D.x)
                        der isa NoDeriv && return D(arg)
                        push!(summed_args, der[idx])
                        continue
                    elseif der isa NoDeriv
                        push!(summed_args, Differential(a)(arg) * executediff(D, a))
                    else
                        push!(summed_args, der[idx] * executediff(D, a))
                    end
                end
                return SymbolicUtils.add_worker(VartypeT, summed_args)
            elseif f === ifelse
                inner_args = arguments(arg)
                dtrue = executediff(D, inner_args[2]; throw_no_derivative)
                dfalse = executediff(D, inner_args[3]; throw_no_derivative)
                args = SymbolicUtils.ArgsT{VartypeT}((inner_args[1], dtrue, dfalse))
                return BSImpl.Term{VartypeT}(ifelse, args; type = symtype(arg), shape = shape(arg))
            elseif f isa Differential
                # The recursive expand_derivatives was not able to remove
                # a nested Differential. We can attempt to differentiate the
                # inner expression wrt to the outer iv. And leave the
                # unexpandable Differential outside.
                isequal(f.x, D.x) && return D(arg)
                inner_args = arguments(arg)
                innerdiff = executediff(D, inner_args[1]; simplify, throw_no_derivative)
                return @match innerdiff begin
                    BSImpl.Term(; f = finner) && if finner isa Differential end => D(arg)
                    _ => executediff(f, innerdiff; simplify, throw_no_derivative)
                end
            elseif f isa Integral && f.domain.domain isa AbstractInterval
                domain = f.domain.domain
                domainvars = f.domain.variables
                a, b = unwrap.(DomainSets.endpoints(domain))
                summed_args = SymbolicUtils.ArgsT{VartypeT}()
                inner_function = arguments(arg)[1]
                if iscall(a) || isequal(a, D.x)
                    t1 = substitute_in_deriv(inner_function, Dict(domainvars => a))
                    t2 = executediff(D, a; simplify, throw_no_derivative)
                    push!(summed_args, -t1*t2)
                end
                if iscall(b) || isequal(b, D.x)
                    t1 = substitute_in_deriv(inner_function, Dict(domainvars => b))
                    t2 = executediff(D, b; simplify, throw_no_derivative)
                    push!(summed_args, t1*t2)
                end
                inner = executediff(D, inner_function; simplify, throw_no_derivative)
                push!(summed_args, f(inner))
                return SymbolicUtils.add_worker(VartypeT, summed_args)
            elseif f === (^) && SymbolicUtils.isconst(args[2])
                base, exp = args
                prod_args = (exp, (base ^ Const{VartypeT}(exp - 1))::BasicSymbolic{VartypeT}, executediff(D, base; simplify, throw_no_derivative))
                return SymbolicUtils.mul_worker(VartypeT, prod_args)
            else
                inner_args = arguments(arg)
                summed_args = SymbolicUtils.ArgsT{VartypeT}()

                for (i, iarg) in enumerate(inner_args)
                    t2 = executediff(D, iarg; simplify, throw_no_derivative)::SymbolicT
                    _iszero(t2) && continue
                    t = derivative_idx(arg, i)::Union{NoDeriv, SymbolicT}
                    if t isa NoDeriv
                        throw_no_derivative && throw(DerivativeNotDefinedError(arg, i))
                        t = D(arg)
                    end
                    if !_isone(t2)
                        t = t * t2
                    end
                    push!(summed_args, t)
                end
                return SymbolicUtils.add_worker(VartypeT, summed_args)
            end
        end
        BSImpl.AddMul(; coeff, dict, variant) => begin
                if variant == SymbolicUtils.AddMulVariant.ADD
                    inner_args = arguments(arg)
                    summed_args = SymbolicUtils.ArgsT{VartypeT}()
                    for iarg in inner_args
                        t2 = executediff(D, iarg; simplify, throw_no_derivative)
                        _iszero(t2) && continue
                        push!(summed_args, t2)
                    end
                    return SymbolicUtils.add_worker(VartypeT, summed_args)
                else
                    # Do the `add_with_div` trick where we write to `inner_args` to avoid
                    # copying the array but restore its state after.
                    # TODO: This might be possible to do faster by using `_mul_worker!`
                    # directly
                    inner_args = parent(arguments(arg))
                    summed_args = SymbolicUtils.ArgsT{VartypeT}()

                    for (i, iarg) in enumerate(inner_args)
                            t2 = executediff(D, iarg; simplify, throw_no_derivative)
                            _iszero(t2) && continue
                        try
                            inner_args[i] = t2
                            push!(summed_args, SymbolicUtils.mul_worker(VartypeT, inner_args))
                        finally
                            inner_args[i] = iarg
                        end
                    end
                    return SymbolicUtils.add_worker(VartypeT, summed_args)
                end
        end
        BSImpl.Div(; num, den) => begin
            dnum = executediff(D, num; simplify, throw_no_derivative)
            dden = executediff(D, den; simplify, throw_no_derivative)
            return (dnum / den - num * dden / den^2)
            newnum = den * dnum - num * dden
            newden = dden ^ 2
            return SymbolicUtils.Div{VartypeT}(newnum, newden, false; type = symtype(arg), shape = shape(arg))
        end
        BSImpl.ArrayOp(; output_idx, expr, reduce, term, ranges) => begin
            if term !== nothing
                term = executediff(D, term; simplify, throw_no_derivatives)
            end
            expr = executediff(D, expr; simplify, throw_no_derivatives)
            return BSImpl.ArrayOp{VartypeT}(output_idx, expr, reduce, term, ranges; type = symtype(arg), shape = shape(arg))
        end
    end
end

"""
$(SIGNATURES)

Expands derivatives within a symbolic expression `O`.

This function recursively traverses a symbolic expression, applying the chain rule
and other derivative rules to expand any derivatives it encounters.

# Arguments
- `O::BasicSymbolic`: The symbolic expression to expand.
- `simplify::Bool=false`: Whether to simplify the resulting expression using
    [`SymbolicUtils.simplify`](@ref).

# Keyword Arguments
- `throw_no_derivative=false`: Whether to throw if a function with unknown
   derivative is encountered.

# Examples
```jldoctest
julia> @variables x y z k;

julia> f = k*(abs(x-y)/y-z)^2
k*((abs(x - y) / y - z)^2)

julia> Dx = Differential(x) # Differentiate wrt x
(::Differential) (generic function with 2 methods)

julia> dfx = expand_derivatives(Dx(f))
(k*((2abs(x - y)) / y - 2z)*ifelse(signbit(x - y), -1, 1)) / y
```
"""
function expand_derivatives(O::BasicSymbolic, simplify=false; throw_no_derivative=false)
    @match O begin
        BSImpl.Term(; f, args) && if f isa Differential end => begin
            arg = expand_derivatives(args[1], false; throw_no_derivative)
            return executediff(f, arg; simplify, throw_no_derivative)
        end
        BSImpl.Term(; f, args) && if f isa Integral end => begin
            return f(expand_derivatives(args[1]; throw_no_derivative))
        end
        if !hasderiv(O) end => return O
        _ => begin
            newargs = SArgsT()
            for arg in arguments(O)
                push!(newargs, expand_derivatives(arg, false; throw_no_derivative))
            end
            O1 = maketerm(SymbolicT, operation(O), newargs, nothing; type = symtype(O))
            return simplify ? SymbolicUtils.simplify(O1) : O1
        end
    end
end
function expand_derivatives(n::Num, simplify=false; kwargs...)
    Num(expand_derivatives(value(n), simplify; kwargs...))
end
function expand_derivatives(n::Complex{Num}, simplify=false; kwargs...)
    re = expand_derivatives(real(n), simplify; kwargs...)
    img = expand_derivatives(imag(n), simplify; kwargs...)
    Complex{Num}(Num(re), Num(img))
end
expand_derivatives(x, simplify=false; kwargs...) = x

# Don't specialize on the function here
"""
$(SIGNATURES)

Calculate the derivative of the op `O` with respect to its argument with index
`idx`.

# Examples

```jldoctest label1
julia> using Symbolics

julia> @variables x y;

julia> Symbolics.derivative_idx(Symbolics.value(sin(x)), 1)
cos(x)
```

Note that the function does not recurse into the operation's arguments, i.e., the
chain rule is not applied:

```jldoctest label1
julia> myop = Symbolics.value(sin(x) * y^2)
sin(x)*(y^2)

julia> typeof(Symbolics.operation(myop))  # Op is multiplication function
typeof(*)

julia> Symbolics.derivative_idx(myop, 1)  # wrt. sin(x)
y^2

julia> Symbolics.derivative_idx(myop, 2)  # wrt. y^2
sin(x)
```
"""
derivative_idx(O::Any, ::Any) = COMMON_ZERO
function derivative_idx(O::BasicSymbolic, idx)
    iscall(O) || return COMMON_ZERO
    res = derivative(operation(O), (arguments(O)...,), Val(idx))
    if res isa NoDeriv
        return res
    else
        return Const{VartypeT}(res)
    end
end

# Indicate that no derivative is defined.
struct NoDeriv
end

"""
    Symbolics.derivative(::typeof(f), args::NTuple{N, Any}, ::Val{i})

Return the derivative of `f(args...)` with respect to `args[i]`. `N` should be the number
of arguments that `f` takes and `i` is the argument with respect to which the derivative
is taken. The result can be a numeric value (if the derivative is constant) or a symbolic
expression. This function is useful for defining derivatives of custom functions registered
via `@register_symbolic`, to be used when calling `expand_derivatives`.
"""
derivative(f, args, v) = NoDeriv()

# Pre-defined derivatives
import DiffRules
for (modu, fun, arity) ∈ DiffRules.diffrules(; filter_modules=(:Base, :SpecialFunctions, :NaNMath))
    fun in [:*, :+, :abs, :mod, :rem, :max, :min] && continue # special
    for i ∈ 1:arity

        expr = if arity == 1
            DiffRules.diffrule(modu, fun, :(args[1]))
        else
            DiffRules.diffrule(modu, fun, ntuple(k->:(args[$k]), arity)...)[i]
        end
        @eval derivative(::typeof($modu.$fun), args::NTuple{$arity,Any}, ::Val{$i}) = $expr
    end
end

derivative(::typeof(+), args::NTuple{N,Any}, ::Val) where {N} = 1
derivative(::typeof(*), args::NTuple{N,Any}, ::Val{i}) where {N,i} = *(deleteat!(collect(args), i)...)
derivative(::typeof(one), args::Tuple{<:Any}, ::Val) = 0

derivative(f::Function, x::Union{Num, <:BasicSymbolic}) = derivative(f(x), x)
derivative(::Function, x::Any) = TypeError(:derivative, "2nd argument", Union{Num, <:BasicSymbolic}, x) |> throw

function count_order(x)
    @assert !(x isa Symbol) "The variable $x must have an order of differentiation that is greater or equal to 1!"
    n = 1
    while !(x.args[1] isa Symbol)
        n = n+1
        x = x.args[1]
    end
    n, x.args[1]
end

_repeat_apply(f, n) = n == 1 ? f : ComposedFunction{Any,Any}(f, _repeat_apply(f, n-1))
function _differential_macro(x)
    ex = Expr(:block)
    push!(ex.args,  :(Base.depwarn("`@derivatives D'''~x` is deprecated. Use `Differential(x)^3` instead.", Symbol("@derivatives"))))
    lhss = Symbol[]
    x = x isa Tuple && first(x).head == :tuple ? first(x).args : x # tuple handling
    x = flatten_expr!(x)
    for di in x
        @assert di isa Expr && di.args[1] == :~ "@derivatives expects a form that looks like `@derivatives D''~t E'~t` or `@derivatives (D''~t), (E'~t)`"
        lhs = di.args[2]
        rhs = di.args[3]
        order, lhs = count_order(lhs)
        push!(lhss, lhs)
        expr = :($lhs = Differential($rhs, $order))
        push!(ex.args,  expr)
    end
    push!(ex.args, Expr(:tuple, lhss...))
    ex
end

"""
$(SIGNATURES)

Define one or more differentials.

# Examples

```jldoctest
julia> using Symbolics

julia> @variables x y z;

julia> Dx = Differential(x); Dy = Differential(y);  # Create differentials wrt. x and y

julia> Dx(z)  # Differentiate z wrt. x
Differential(x)(z)

julia> Dy(z)  # Differentiate z wrt. y
Differential(y)(z)
```
"""
macro derivatives(x...)
    esc(_differential_macro(x))
end


### Jacobians & Hessians

"""
$(SIGNATURES)

A helper function for computing the derivative of the expression `O` with respect to
`var`.

# Keyword Arguments

- `simplify=false`: The simplify argument of `expand_derivatives`.

All other keyword arguments are forwarded to `expand_derivatives`.
"""
function derivative(O, var; simplify=false, kwargs...)
    if O isa AbstractArray
        Num[Num(expand_derivatives(Differential(var)(unwrap(o)), simplify; kwargs...)) for o in O]
    else
        Num(expand_derivatives(Differential(var)(unwrap(O)), simplify; kwargs...))
    end
end

"""
$(SIGNATURES)

A helper function for computing the gradient of the expression `O` with respect to
an array of variable expressions.

# Keyword Arguments

- `simplify=false`: The simplify argument of `expand_derivatives`.

All other keyword arguments are forwarded to `expand_derivatives`.
"""
function gradient(O, vars::AbstractVector; simplify=false, kwargs...)
    Num[Num(expand_derivatives(Differential(vars[vi])(unwrap(O)),simplify; kwargs...)) for vi in eachindex(vars)]
end

"""
$(SIGNATURES)

A helper function for computing the Jacobian of an array of expressions with respect to
an array of variable expressions.

# Keyword Arguments

- `simplify=false`: The simplify argument of `expand_derivatives`.
- `scalarize=true`: Whether to scalarize `ops` and `vars` before computing the jacobian.

All other keyword arguments are forwarded to `expand_derivatives`.
"""
function jacobian(ops::AbstractVector, vars::AbstractVector; simplify=false, scalarize=true, kwargs...)
    if scalarize
        ops = Symbolics.scalarize(ops)
        vars = Symbolics.scalarize(vars)
    end
    Num[Num(expand_derivatives(Differential(unwrap(v))(unwrap(O)),simplify; kwargs...)) for O in ops, v in vars]
end

function jacobian(ops, vars; simplify=false, kwargs...)
    ops = vec(scalarize(ops))
    vars = vec(scalarize(vars)) # Suboptimal, but prevents wrong results on Arr for now. Arr resulting from a symbolic function will fail on this due to unknown size.
    jacobian(ops, vars; simplify=simplify, scalarize=false, kwargs...)
end

function faster_maybe_scalarize!(arg::Vector)
    for (i, x) in enumerate(arg)
        arg[i] = scalarize(x)
    end
    return arg
end

faster_maybe_scalarize!(arg) = scalarize(arg)

"""
$(SIGNATURES)

A helper function for computing the sparse Jacobian of an array of expressions with respect to
an array of variable expressions.

# Keyword Arguments

- `simplify=false`: The simplify argument of `expand_derivatives`.

All other keyword arguments are forwarded to `expand_derivatives`.
"""
function sparsejacobian(ops::AbstractVector, vars::AbstractVector; simplify::Bool=false, kwargs...)
    ops = faster_maybe_scalarize!(ops)
    vars = faster_maybe_scalarize!(vars)
    sp = jacobian_sparsity(ops, vars)
    I,J,_ = findnz(sp)

    exprs = sparsejacobian_vals(ops, vars, I, J; simplify=simplify, kwargs...)

    sparse(I, J, exprs, length(ops), length(vars))
end

"""
$(SIGNATURES)

A helper function for computing the values of the sparse Jacobian of an array of expressions with respect to
an array of variable expressions given the sparsity structure.

# Keyword Arguments

- `simplify=false`: The simplify argument of `expand_derivatives`.

All other keyword arguments are forwarded to `expand_derivatives`.
"""
function sparsejacobian_vals(ops::AbstractVector, vars::AbstractVector, I::AbstractVector, J::AbstractVector; simplify::Bool=false, kwargs...)
    ops = faster_maybe_scalarize!(ops)
    vars = faster_maybe_scalarize!(vars)

    exprs = Num[]
    sizehint!(exprs, length(I))

    for (i,j) in zip(I, J)
        push!(exprs, Num(expand_derivatives(Differential(vars[j])(ops[i]), simplify; kwargs...)))
    end
    exprs
end

"""
$(TYPEDSIGNATURES)

Return the sparsity pattern of the Jacobian of an array of expressions with respect to
an array of variable expressions.

# Arguments
- `exprs`: an array of symbolic expressions.
- `vars`: an array of symbolic variables.

# Examples
```jldoctest
julia> using Symbolics

julia> vars = @variables x₁ x₂;

julia> exprs = [2x₁, 3x₂, 4x₁ * x₂];

julia> Symbolics.jacobian_sparsity(exprs, vars)
3×2 SparseArrays.SparseMatrixCSC{Bool, Int64} with 4 stored entries:
 1  ⋅
 ⋅  1
 1  1
```
"""
function jacobian_sparsity(exprs::AbstractArray, vars::AbstractArray)
    if any(iswrapped, exprs)
        du = map(value, exprs)
    else
        du = exprs
    end
    if any(iswrapped, vars)
        u = map(value, vars)
    else
        u = vars
    end
    dict = Dict(zip(u, 1:length(u)))

    i = Ref(1)
    I = Int[]
    J = Int[]
    sizehint!(I, 2length(exprs))
    sizehint!(J, 2length(vars))


    # This rewriter notes down which u's appear in a
    # given du (whose index is stored in the `i` Ref)

    function r(x)
        if iscall(x)
            for y in arguments(x)
                r(y)
            end
        end
        j = get(dict, x, -1)
        if j != -1
            push!(I, i[])
            push!(J, j)
        end
    end

    for ii = 1:length(du)
        i[] = ii
        r(du[ii])
    end

    sparse(I, J, true, length(du), length(u))
end
"""
$(TYPEDSIGNATURES)

Return the sparsity pattern of the Jacobian of the mutating function `f!`.

# Arguments
- `f!`: an in-place function `f!(output, input, args...; kwargs...)`.
- `output`: output array.
- `input`: input array.

The [eltype](https://docs.julialang.org/en/v1/base/collections/#Base.eltype)
of `output` and `input` can be either symbolic or
[primitive](https://docs.julialang.org/en/v1/manual/types/#Primitive-Types).

# Examples
```jldoctest
julia> using Symbolics

julia> f!(y, x) = y .= [x[2], 2x[1], 3x[1] * x[2]];

julia> output = Vector{Float64}(undef, 3);

julia> input = Vector{Float64}(undef, 2);

julia> Symbolics.jacobian_sparsity(f!, output, input)
3×2 SparseArrays.SparseMatrixCSC{Bool, Int64} with 4 stored entries:
 ⋅  1
 1  ⋅
 1  1
```
"""
function jacobian_sparsity(f!, output::AbstractArray, input::AbstractArray,
                           args...; kwargs...)
    exprs = similar(output, Num)
    fill!(exprs, false)
    vars = ArrayInterface.restructure(input, map(variable, eachindex(input)))
    f!(exprs, vars, args...; kwargs...)
    jacobian_sparsity(exprs, vars)
end

"""
    exprs_occur_in(exprs::Vector, expr)

Return an array of booleans `finds` where `finds[i]` is true if `exprs[i]` occurs in `expr`
false otherwise.
"""
function exprs_occur_in(exprs, expr)
    vec(jacobian_sparsity([expr], exprs))
end

"""
$(SIGNATURES)

A helper function for computing the Hessian of the expression `O` with respect to
an array of variable expressions.

# Keyword Arguments

- `simplify=false`: The simplify argument of `expand_derivatives`.

All other keyword arguments are forwarded to `expand_derivatives`.
"""
function hessian(O, vars::AbstractVector; simplify=false, kwargs...)
    vars = map(value, vars)
    first_derivs = map(value, vec(jacobian([values(O)], vars; simplify=simplify, kwargs...)))
    n = length(vars)
    H = Array{Num, 2}(undef,(n, n))
    fill!(H, 0)
    for i=1:n
        for j=1:i
            H[j, i] = H[i, j] = expand_derivatives(Differential(vars[i])(first_derivs[j]), simplify; kwargs...)
        end
    end
    H
end

hessian(O, vars::Arr; kwargs...) = hessian(O, collect(vars); kwargs...) 

isidx(x) = unwrap_const(x) isa TermCombination

basic_mkterm(t, g, args, m) = metadata(Term{VartypeT}(g, args; type = Any), m)

const _scalar = one(TermCombination)

const linearity_rules = (
      (@rule +(~~xs) => reduce(+, filter(isidx, map(unwrap_const, ~~xs)), init=_scalar)),
      (@rule *(~~xs) => reduce(*, filter(isidx, map(unwrap_const, ~~xs)), init=_scalar)),

      (@rule (~f)(~x) => isidx(~x) ? combine_terms_1(linearity_1(~f), ~x) : _scalar),
      (@rule (^)(~x::isidx, ~y) => ~y isa Number && isone(~y) ? ~x : (~x) * (~x)),
      (@rule (~f)(~x, ~y) => combine_terms_2(linearity_2(~f), isidx(~x) ? ~x : _scalar, isidx(~y) ? ~y : _scalar)),

      (@rule ~x::issym => 0),

      # `ifelse(cond, x, y)` can be written as cond * x + (1 - cond) * y
      # where condition `cond` is considered constant in differentiation
      (@rule ifelse(~cond, ~x, ~y) => (isidx(~x) ? ~x : _scalar) + (isidx(~y) ? ~y : _scalar)),

      # Fallback: Unknown functions with arbitrary number of arguments have non-zero partial derivatives
      # Functions with 1 and 2 arguments are already handled above
      (@rule (~f)(~~xs) => reduce(+, filter(isidx, map(unwrap_const, ~~xs)); init=_scalar)^2),
)
const linearity_rules_affine = (
      (@rule +(~~xs) => reduce(+, filter(isidx, map(unwrap_const, ~~xs)), init=_scalar)),
      (@rule *(~~xs) => reduce(*, filter(isidx, map(unwrap_const, ~~xs)), init=_scalar)),

      (@rule (~f)(~x) => isidx(~x) ? combine_terms_1(linearity_1(~f), ~x) : _scalar),
      (@rule (^)(~x::isidx, ~y) => ~y isa Number && isone(~y) ? unwrap_const(~x) : unwrap_const(~x) * unwrap_const(~x)),
      (@rule (~f)(~x, ~y) => combine_terms_2(linearity_2(~f), isidx(~x) ? unwrap_const(~x) : _scalar, isidx(~y) ? unwrap_const(~y) : _scalar)),

      (@rule ~x::issym => 0),
      # if the condition is dependent on the variable, do not consider this as affine
      (@rule ifelse(~cond::isidx, ~x, ~y) => (~cond)^2),
      # `ifelse(cond, x, y)` can be written as cond * x + (1 - cond) * y
      # where condition `cond` is considered constant in differentiation
      (@rule ifelse(~cond::(!isidx), ~x, ~y) => (isidx(~x) ? unwrap_const(~x) : _scalar) + (isidx(~y) ? unwrap_const(~y) : _scalar)),
      # Fallback: Unknown functions with arbitrary number of arguments have non-zero partial derivatives
      # Functions with 1 and 2 arguments are already handled above
      (@rule (~f)(~~xs) => reduce(+, filter(isidx, map(unwrap_const, ~~xs)); init=_scalar)^2),
)
const linearity_propagator = Fixpoint(Postwalk(Chain(linearity_rules); maketerm=basic_mkterm))
const affine_linearity_propagator = Fixpoint(Postwalk(Chain(linearity_rules_affine); maketerm=basic_mkterm))

"""
$(TYPEDSIGNATURES)

Return the sparsity pattern of the Hessian of an expression with respect to
an array of variable expressions.

# Arguments
- `expr`: a symbolic expression.
- `vars`: a vector of symbolic variables.

# Examples
```jldoctest
julia> using Symbolics

julia> vars = @variables x₁ x₂;

julia> expr = 3x₁^2 + 4x₁ * x₂;

julia> Symbolics.hessian_sparsity(expr, vars)
2×2 SparseArrays.SparseMatrixCSC{Bool, Int64} with 3 stored entries:
 1  1
 1  ⋅
```
"""
function hessian_sparsity(expr, vars::AbstractVector; full::Bool=true, linearity_propagator = linearity_propagator)
    @assert !(expr isa AbstractArray)
    expr = value(expr)
    u = map(value, vars)
    dict = Dict(ui => TermCombination(Set([Dict(i=>1)])) for (i, ui) in enumerate(u))
    f = Rewriters.Prewalk(x-> get(dict, x, x); maketerm=basic_mkterm)(expr)
    lp = unwrap_const(linearity_propagator(f))
    S = _sparse(lp, length(u))
    S = full ? S : tril(S)
end

"""
$(TYPEDSIGNATURES)

Return the sparsity pattern of the Hessian of the given function `f`.

# Arguments
- `f`: an out-of-place function `f(input, args...; kwargs...)`.
- `input`: a vector of input values whose [eltype](https://docs.julialang.org/en/v1/base/collections/#Base.eltype) can be either symbolic or [primitive](https://docs.julialang.org/en/v1/manual/types/#Primitive-Types).

# Examples
```jldoctest
julia> using Symbolics

julia> f(x) = 4x[1] * x[2] - 5x[2]^2;

julia> input = Vector{Float64}(undef, 2);

julia> Symbolics.hessian_sparsity(f, input)
2×2 SparseArrays.SparseMatrixCSC{Bool, Int64} with 3 stored entries:
 ⋅  1
 1  1
```
"""
function hessian_sparsity(f::Function, input::AbstractVector, args...; full::Bool=true, kwargs...)
    vars = ArrayInterface.restructure(input, map(variable, eachindex(input)))
    expr = f(vars, args...; kwargs...)
    hessian_sparsity(expr, vars, full=full)
end

"""
$(SIGNATURES)

Check if an expression is affine with respect to a list of variable expressions.
"""
function isaffine(ex, u)
    isempty(hessian_sparsity(ex, u; linearity_propagator = affine_linearity_propagator).nzval)
end

"""
$(SIGNATURES)

Check if an expression is linear with respect to a list of variable expressions.
"""
function islinear(ex, u)
    isaffine(ex, u) && iszero(Num(substitute(ex, Dict(u .=> 0))))
end

"""
$(SIGNATURES)

A helper function for computing the sparse Hessian of an expression with respect to
an array of variable expressions.

# Keyword Arguments

- `simplify=false`: The simplify argument of `expand_derivatives`.
- `full=false`: Whether to construct the full hessian by also including entries in
  the upper-triangular half of the matrix.

All other keyword arguments are forwarded to `expand_derivatives`.
"""
function sparsehessian(op, vars::AbstractVector; simplify::Bool=false, full::Bool=true, kwargs...)
    op = value(op)
    vars = map(value, vars)
    S = hessian_sparsity(op, vars, full=full)
    I, J, _ = findnz(S)

    exprs = sparsehessian_vals(op, vars, I, J; simplify=simplify, kwargs...)

    H = sparse(I, J, exprs, length(vars), length(vars))

    if full
        for (i, j) in zip(I, J)
            j > i && (H[i, j] = H[j, i])
        end
    end
    return H
end

"""
$(SIGNATURES)

A helper function for computing the values of the sparse Hessian of an expression with respect to
an array of variable expressions given the sparsity structure.

# Keyword Arguments

- `simplify=false`: The simplify argument of `expand_derivatives`.

All other keyword arguments are forwarded to `expand_derivatives`.
"""
function sparsehessian_vals(op, vars::AbstractVector, I::AbstractVector, J::AbstractVector; simplify::Bool=false, kwargs...)
    vars = Symbolics.scalarize(vars)

    exprs = Array{Num}(undef, length(I))
    fill!(exprs, 0)

    prev_j = 0
    d = nothing
    for (k, (i, j)) in enumerate(zip(I, J))
        j > i && continue
        if j != prev_j
            d = expand_derivatives(Differential(vars[j])(op), false; kwargs...)
        end
        expr = expand_derivatives(Differential(vars[i])(d), simplify; kwargs...)
        exprs[k] = expr
        prev_j = j
    end
    exprs
end
