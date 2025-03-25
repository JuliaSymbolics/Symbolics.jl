abstract type Operator end
propagate_shape(::Operator, x) = axes(x)

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
    x
    Differential(x) = new(value(x))
end
function (D::Differential)(x)
    x = unwrap(x)
    if isarraysymbolic(x)
        array_term(D, x)
    else
        term(D, x)
    end
end
(D::Differential)(x::Union{Num, Arr}) = wrap(D(unwrap(x)))
(D::Differential)(x::Complex{Num}) = wrap(ComplexTerm{Real}(D(unwrap(real(x))), D(unwrap(imag(x)))))
SymbolicUtils.promote_symtype(::Differential, T) = T
SymbolicUtils.isbinop(f::Differential) = false

is_derivative(x) = iscall(x) ? operation(x) isa Differential : false

Base.:*(D1, D2::Differential) = D1 ∘ D2
Base.:*(D1::Differential, D2) = D1 ∘ D2
Base.:*(D1::Differential, D2::Differential) = D1 ∘ D2
Base.:^(D::Differential, n::Integer) = iszero(n) ? identity : _repeat_apply(D, n)

Base.show(io::IO, D::Differential) = print(io, "Differential(", D.x, ")")
Base.nameof(D::Differential) = :Differential

Base.:(==)(D1::Differential, D2::Differential) = isequal(D1.x, D2.x)
Base.hash(D::Differential, u::UInt) = hash(D.x, xor(u, 0xdddddddddddddddd))

_isfalse(occ::Bool) = occ === false
_isfalse(occ::Symbolic) = iscall(occ) && _isfalse(operation(occ))


SymbolicUtils.@cache limit = 500_000 function occursin_info(x::BasicSymbolic, expr::Any, fail::Bool = true)::Bool
    _occursin_info(x, expr, fail)
end

function _occursin_info(x, expr, fail = true)
    if symtype(expr) <: AbstractArray
        if fail
            error("Differentiation with array expressions is not yet supported")
        else
            return occursin(x, expr)
        end
    end

    # Allow scalarized expressions
    function is_scalar_indexed(ex)
        (iscall(ex) && operation(ex) == getindex && !(symtype(ex) <: AbstractArray)) ||
        (iscall(ex) && (issym(operation(ex)) || iscall(operation(ex))) &&
         is_scalar_indexed(operation(ex)))
    end

    # x[1] == x[1] but not x[2]
    if is_scalar_indexed(x) && is_scalar_indexed(expr) &&
        isequal(first(arguments(x)), first(arguments(expr)))
        return isequal(operation(x), operation(expr)) &&
               isequal(arguments(x), arguments(expr))
    end

    if is_scalar_indexed(x) && is_scalar_indexed(expr) &&
        !occursin(first(arguments(x)), first(arguments(expr)))
        return false
    end

    if is_scalar_indexed(expr) && !is_scalar_indexed(x) && !occursin(x, expr)
        return false
    end

    !iscall(expr) && return isequal(x, expr)
    if isequal(x, expr)
        true
    else
        any(a -> occursin_info(x, a, operation(expr) !== getindex), arguments(expr))
    end
end

function _occursin_info(x, expr::Sym, fail)
    if symtype(expr) <: AbstractArray && fail
            error("Differentiation of expressions involving arrays and array variables is not yet supported.")
    end
    isequal(x, expr)
end

"""
    hasderiv(O)

Returns true if the expression or equation `O` contains [`Differential`](@ref) terms.
"""
hasderiv(O) = recursive_hasoperator(Differential, O)


_recursive_hasoperator(op, eq::Equation) = recursive_hasoperator(op, eq.lhs) || recursive_hasoperator(op, eq.rhs)
_recursive_hasoperator(op) = Base.Fix1(recursive_hasoperator, op) # curry version
_recursive_hasoperator(::Type{T}, ::T) where T = true


"""
    recursive_hasoperator(op, O)

An internal function that contains the logic for [`hasderiv`](@ref) and [`hasdiff`](@ref).
Return true if `O` contains a term with `Operator` `op`.
"""
SymbolicUtils.@cache function recursive_hasoperator(op::Any, O::Any)::Bool
    _recursive_hasoperator(op, O)
end

function _recursive_hasoperator(op, O)
    iscall(O) || return false
    if operation(O) isa op
        return true
    else
        if isadd(O) || ismul(O)
            any(_recursive_hasoperator(op), keys(O.dict))
        elseif ispow(O)
            _recursive_hasoperator(op)(O.base) || _recursive_hasoperator(op)(O.exp)
        elseif isdiv(O)
            _recursive_hasoperator(op)(O.num) || _recursive_hasoperator(op)(O.den)
        else
            any(_recursive_hasoperator(op), arguments(O))
        end
    end
end

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

"""
    executediff(D, arg, simplify=false; occurrences=nothing)

Apply the passed Differential D on the passed argument.

This function differs to `expand_derivatives` in that in only expands the
passed differential and not any other Differentials it encounters.

# Arguments
- `D::Differential`: The differential to apply
- `arg::Symbolic`: The symbolic expression to apply the differential on.
- `simplify::Bool=false`: Whether to simplify the resulting expression using
    [`SymbolicUtils.simplify`](@ref).
- `occurrences=nothing`: Information about the occurrences of the independent
    variable in the argument of the derivative. This is used internally for
    optimization purposes.
- `throw_no_derivative=false`: Whether to throw if a function with unknown
    derivative is encountered.
"""
function executediff(D, arg, simplify=false; throw_no_derivative=false)
    isequal(arg, D.x) && return 1
    occursin_info(D.x, arg) || return 0

    if !iscall(arg)
        return D(arg) # Cannot expand
    elseif (op = operation(arg); issym(op))
        inner_args = arguments(arg)
        if any(isequal(D.x), inner_args)
            return D(arg) # base case if any argument is directly equal to the i.v.
        else
            return sum(inner_args, init=0) do a
                return executediff(Differential(a), arg; throw_no_derivative) *
                executediff(D, a; throw_no_derivative)
            end
        end
    elseif op === getindex
        inner_args = arguments(arguments(arg)[1])
        c = 0
        for a in inner_args
            if isequal(a, D.x)
                return D(arg)
            else
                c += Differential(a)(arg) * D(a)
            end
        end
        return expand_derivatives(c)
    elseif op === ifelse
        args = arguments(arg)
        O = op(args[1], 
            executediff(D, args[2], simplify; throw_no_derivative),
            executediff(D, args[3], simplify; throw_no_derivative))
        return O
    elseif isa(op, Differential)
        # The recursive expand_derivatives was not able to remove
        # a nested Differential. We can attempt to differentiate the
        # inner expression wrt to the outer iv. And leave the
        # unexpandable Differential outside.
        if isequal(op.x, D.x)
            return D(arg)
        else
            inner = executediff(D, arguments(arg)[1], false; throw_no_derivative)
            # if the inner expression is not expandable either, return
            if iscall(inner) && operation(inner) isa Differential
                return D(arg)
            else
                # otherwise give the nested Differential another try
                return executediff(op, inner, simplify; throw_no_derivative)
            end
        end
    elseif isa(op, Integral)
        if isa(op.domain.domain, AbstractInterval)
            domain = op.domain.domain
            a, b = DomainSets.endpoints(domain)
            c = 0
            inner_function = arguments(arg)[1]
            if iscall(value(a))
                t1 = SymbolicUtils.substitute(inner_function, Dict(op.domain.variables => value(a)))
                t2 = D(a)
                c -= t1*t2
            end
            if iscall(value(b))
                t1 = SymbolicUtils.substitute(inner_function, Dict(op.domain.variables => value(b)))
                t2 = D(b)
                c += t1*t2
            end
            inner = executediff(D, arguments(arg)[1]; throw_no_derivative)
            c += op(inner)
            return value(c)
        end
    end

    inner_args = arguments(arg)
    l = length(inner_args)
    exprs = []
    c = 0

    for i in 1:l
        t2 = executediff(D, inner_args[i],false; throw_no_derivative)

        x = if _iszero(t2)
            t2
        elseif _isone(t2)
            d = derivative_idx(arg, i)
            if d isa NoDeriv
                throw_no_derivative && throw(DerivativeNotDefinedError(arg, i))
                D(arg)
            else
                d
            end
        else
            t1 = derivative_idx(arg, i)
            t1 = if t1 isa NoDeriv
                throw_no_derivative && throw(DerivativeNotDefinedError(arg, i))
                D(arg)
            else
                t1
            end
            t1 * t2
        end

        if _iszero(x)
            continue
        elseif x isa Symbolic
            push!(exprs, x)
        else
            c += x
        end
    end

    if isempty(exprs)
        return c
    elseif length(exprs) == 1
        term = (simplify ? SymbolicUtils.simplify(exprs[1]) : exprs[1])
        return _iszero(c) ? term : c + term
    else
        x = +((!_iszero(c) ? vcat(c, exprs) : exprs)...)
        return simplify ? SymbolicUtils.simplify(x) : x
    end
end

"""
$(SIGNATURES)

Expands derivatives within a symbolic expression `O`.

This function recursively traverses a symbolic expression, applying the chain rule
and other derivative rules to expand any derivatives it encounters.

# Arguments
- `O::Symbolic`: The symbolic expression to expand.
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
function expand_derivatives(O::Symbolic, simplify=false; throw_no_derivative=false)
    if iscall(O) && isa(operation(O), Differential)
        arg = only(arguments(O))
        arg = expand_derivatives(arg, false; throw_no_derivative)
        return executediff(operation(O), arg, simplify; throw_no_derivative)
    elseif iscall(O) && isa(operation(O), Integral)
        return operation(O)(expand_derivatives(arguments(O)[1]; throw_no_derivative))
    elseif !hasderiv(O)
        return O
    else
        args = map(a->expand_derivatives(a, false; throw_no_derivative), arguments(O))
        O1 = operation(O)(args...)
        return simplify ? SymbolicUtils.simplify(O1) : O1
    end
end
function expand_derivatives(n::Num, simplify=false; kwargs...)
    wrap(expand_derivatives(value(n), simplify; kwargs...))
end
function expand_derivatives(n::Complex{Num}, simplify=false; kwargs...)
    wrap(ComplexTerm{Real}(expand_derivatives(real(n), simplify; kwargs...),
                           expand_derivatives(imag(n), simplify; kwargs...)))
end
expand_derivatives(x, simplify=false; kwargs...) = x

_iszero(x) = false
_isone(x) = false

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
derivative_idx(O::Any, ::Any) = 0
function derivative_idx(O::Symbolic, idx)
    iscall(O) ? derivative(operation(O), (arguments(O)...,), Val(idx)) : 0
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

derivative(f::Function, x::Num) = derivative(f(x), x)
derivative(::Function, x::Any) = TypeError(:derivative, "2nd argument", Num, typeof(x)) |> throw

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
    push!(ex.args,  :(Base.depwarn("`@derivatives D'''~x` is deprecated. Use `Differential(x)^3` instead.", Symbol("@derivatives"), force=true)))
    lhss = Symbol[]
    x = x isa Tuple && first(x).head == :tuple ? first(x).args : x # tuple handling
    x = flatten_expr!(x)
    for di in x
        @assert di isa Expr && di.args[1] == :~ "@derivatives expects a form that looks like `@derivatives D''~t E'~t` or `@derivatives (D''~t), (E'~t)`"
        lhs = di.args[2]
        rhs = di.args[3]
        order, lhs = count_order(lhs)
        push!(lhss, lhs)
        expr = :($lhs = $_repeat_apply(Differential($value($rhs)), $order))
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
        Num[Num(expand_derivatives(Differential(var)(value(o)), simplify; kwargs...)) for o in O]
    else
        Num(expand_derivatives(Differential(var)(value(O)), simplify; kwargs...))
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
    Num[Num(expand_derivatives(Differential(v)(value(O)),simplify; kwargs...)) for v in vars]
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
    Num[Num(expand_derivatives(Differential(value(v))(value(O)),simplify; kwargs...)) for O in ops, v in vars]
end

function jacobian(ops, vars; simplify=false, kwargs...)
    ops = vec(scalarize(ops))
    vars = vec(scalarize(vars)) # Suboptimal, but prevents wrong results on Arr for now. Arr resulting from a symbolic function will fail on this due to unknown size.
    jacobian(ops, vars; simplify=simplify, scalarize=false, kwargs...)
end

"""
$(SIGNATURES)

A helper function for computing the sparse Jacobian of an array of expressions with respect to
an array of variable expressions.

# Keyword Arguments

- `simplify=false`: The simplify argument of `expand_derivatives`.

All other keyword arguments are forwarded to `expand_derivatives`.
"""
function sparsejacobian(ops::AbstractVector, vars::AbstractVector; simplify::Bool=false, kwargs...)
    ops = Symbolics.scalarize(ops)
    vars = Symbolics.scalarize(vars)
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
    ops = Symbolics.scalarize(ops)
    vars = Symbolics.scalarize(vars)

    exprs = Num[]

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
    du = map(value, exprs)
    u = map(value, vars)
    dict = Dict(zip(u, 1:length(u)))

    i = Ref(1)
    I = Int[]
    J = Int[]


    # This rewriter notes down which u's appear in a
    # given du (whose index is stored in the `i` Ref)

    r = @rule ~x::(x->haskey(dict, x)) => begin
        push!(I, i[])
        push!(J, dict[~x])
        nothing
    end

    r =  Rewriters.Postwalk(r)

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
function jacobian_sparsity(f!::Function, output::AbstractArray, input::AbstractArray,
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

isidx(x) = x isa TermCombination

basic_mkterm(t, g, args, m) = metadata(Term{Any}(g, args), m)

let
    # we do this in a let block so that Revise works on the list of rules

    _scalar = one(TermCombination)

    linearity_rules = [
          @rule +(~~xs) => reduce(+, filter(isidx, ~~xs), init=_scalar)
          @rule *(~~xs) => reduce(*, filter(isidx, ~~xs), init=_scalar)

          @rule (~f)(~x) => isidx(~x) ? combine_terms_1(linearity_1(~f), ~x) : _scalar
          @rule (^)(~x::isidx, ~y) => ~y isa Number && isone(~y) ? ~x : (~x) * (~x)
          @rule (~f)(~x, ~y) => combine_terms_2(linearity_2(~f), isidx(~x) ? ~x : _scalar, isidx(~y) ? ~y : _scalar)

          @rule ~x::issym => 0

          # `ifelse(cond, x, y)` can be written as cond * x + (1 - cond) * y
          # where condition `cond` is considered constant in differentiation
          @rule ifelse(~cond, ~x, ~y) => (isidx(~x) ? ~x : _scalar) + (isidx(~y) ? ~y : _scalar)

          # Fallback: Unknown functions with arbitrary number of arguments have non-zero partial derivatives
          # Functions with 1 and 2 arguments are already handled above
          @rule (~f)(~~xs) => reduce(+, filter(isidx, ~~xs); init=_scalar)^2
    ]
    linearity_propagator = Fixpoint(Postwalk(Chain(linearity_rules); maketerm=basic_mkterm))

    global hessian_sparsity

    @doc """
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
    function hessian_sparsity(expr, vars::AbstractVector; full::Bool=true)
        @assert !(expr isa AbstractArray)
        expr = value(expr)
        u = map(value, vars)
        dict = Dict(ui => TermCombination(Set([Dict(i=>1)])) for (i, ui) in enumerate(u))
        f = Rewriters.Prewalk(x-> get(dict, x, x); maketerm=basic_mkterm)(expr)
        lp = linearity_propagator(f)
        S = _sparse(lp, length(u))
        S = full ? S : tril(S)
    end
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
    isempty(hessian_sparsity(ex, u).nzval)
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

function SymbolicUtils.substitute(op::Differential, dict; kwargs...)
    @set! op.x = substitute(op.x, dict; kwargs...)
end
