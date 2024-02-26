abstract type Operator <: Function end
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

is_derivative(x) = istree(x) ? operation(x) isa Differential : false

Base.:*(D1, D2::Differential) = D1 ∘ D2
Base.:*(D1::Differential, D2) = D1 ∘ D2
Base.:*(D1::Differential, D2::Differential) = D1 ∘ D2
Base.:^(D::Differential, n::Integer) = _repeat_apply(D, n)

Base.show(io::IO, D::Differential) = print(io, "Differential(", D.x, ")")

Base.:(==)(D1::Differential, D2::Differential) = isequal(D1.x, D2.x)
Base.hash(D::Differential, u::UInt) = hash(D.x, xor(u, 0xdddddddddddddddd))

_isfalse(occ::Bool) = occ === false
_isfalse(occ::Symbolic) = istree(occ) && _isfalse(operation(occ))

function occursin_info(x, expr, fail = true)
    if symtype(expr) <: AbstractArray
        if fail
            error("Differentiation with array expressions is not yet supported")
        else
            return occursin(x, expr)
        end
    end

    # Allow scalarized expressions
    function is_scalar_indexed(ex)
        (istree(ex) && operation(ex) == getindex && !(symtype(ex) <: AbstractArray)) ||
        (istree(ex) && (issym(operation(ex)) || istree(operation(ex))) &&
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

    !istree(expr) && return isequal(x, expr)
    if isequal(x, expr)
        true
    else
        args = map(a->occursin_info(x, a, operation(expr) !== getindex), arguments(expr))
        if all(_isfalse, args)
            return false
        end
        Term{Real}(true, args)
    end
end

function occursin_info(x, expr::Sym, fail)
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


recursive_hasoperator(op, eq::Equation) = recursive_hasoperator(op, eq.lhs) || recursive_hasoperator(op, eq.rhs)
recursive_hasoperator(op) = O -> recursive_hasoperator(op, O) # curry version
recursive_hasoperator(::Type{T}, ::T) where T = true


"""
    recursive_hasoperator(op, O)

An internal function that contains the logic for [`hasderiv`](@ref) and [`hasdiff`](@ref).
Return true if `O` contains a term with `Operator` `op`.
"""
function recursive_hasoperator(op, O)
    istree(O) || return false
    if operation(O) isa op
        return true
    else
        if isadd(O) || ismul(O)
            any(recursive_hasoperator(op), keys(O.dict))
        elseif ispow(O)
            recursive_hasoperator(op)(O.base) || recursive_hasoperator(op)(O.exp)
        elseif isdiv(O)
            recursive_hasoperator(op)(O.num) || recursive_hasoperator(op)(O.den)
        else
            any(recursive_hasoperator(op), arguments(O))
        end
    end
end

"""
$(SIGNATURES)

TODO

# Examples
```jldoctest

julia> @variables x y z k;

julia> f=k*(abs(x-y)/y-z)^2
k*((abs(x - y) / y - z)^2)

julia> Dx=Differential(x) # Differentiate wrt x
(::Differential) (generic function with 2 methods)

julia> dfx=expand_derivatives(Dx(f))
(k*((2abs(x - y)) / y - 2z)*IfElse.ifelse(signbit(x - y), -1, 1)) / y
```
"""
function expand_derivatives(O::Symbolic, simplify=false; occurrences=nothing)
    if istree(O) && isa(operation(O), Differential)
        arg = only(arguments(O))
        arg = expand_derivatives(arg, false)

        if occurrences == nothing
            occurrences = occursin_info(operation(O).x, arg)
        end

        _isfalse(occurrences) && return 0
        occurrences isa Bool && return 1 # means it's a `true`

        D = operation(O)

        if !istree(arg)
            return D(arg) # Cannot expand
        elseif (op = operation(arg); issym(op))
            inner_args = arguments(arg)
            if any(isequal(D.x), inner_args)
                return D(arg) # base case if any argument is directly equal to the i.v.
            else
                return sum(inner_args, init=0) do a
                    return expand_derivatives(Differential(a)(arg)) *
                           expand_derivatives(D(a))
                end
            end
        elseif op === (IfElse.ifelse)
            args = arguments(arg)
            O = op(args[1], D(args[2]), D(args[3]))
            return expand_derivatives(O, simplify; occurrences)
        elseif isa(op, Differential)
            # The recursive expand_derivatives was not able to remove
            # a nested Differential. We can attempt to differentiate the
            # inner expression wrt to the outer iv. And leave the
            # unexpandable Differential outside.
            if isequal(op.x, D.x)
                return D(arg)
            else
                inner = expand_derivatives(D(arguments(arg)[1]), false)
                # if the inner expression is not expandable either, return
                if istree(inner) && operation(inner) isa Differential
                    return D(arg)
                else
                    return expand_derivatives(op(inner), simplify)
                end
            end
        elseif isa(op, Integral)
            if isa(op.domain.domain, AbstractInterval)
                domain = op.domain.domain
                a, b = DomainSets.endpoints(domain)
                c = 0
                inner_function = expand_derivatives(arguments(arg)[1])
                if istree(value(a))
                    t1 = SymbolicUtils.substitute(inner_function, Dict(op.domain.variables => value(a)))
                    t2 = D(a)
                    c -= t1*t2
                end
                if istree(value(b))
                    t1 = SymbolicUtils.substitute(inner_function, Dict(op.domain.variables => value(b)))
                    t2 = D(b)
                    c += t1*t2
                end
                inner = expand_derivatives(D(arguments(arg)[1]))
                c += op(inner)
                return value(c)
            end
        end

        inner_args = arguments(arg)
        l = length(inner_args)
        exprs = []
        c = 0

        for i in 1:l
            t2 = expand_derivatives(D(inner_args[i]),false, occurrences=arguments(occurrences)[i])

            x = if _iszero(t2)
                t2
            elseif _isone(t2)
                d = derivative_idx(arg, i)
                d isa NoDeriv ? D(arg) : d
            else
                t1 = derivative_idx(arg, i)
                t1 = t1 isa NoDeriv ? D(arg) : t1
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
    elseif istree(O) && isa(operation(O), Integral)
        return operation(O)(expand_derivatives(arguments(O)[1]))
    elseif !hasderiv(O)
        return O
    else
        args = map(a->expand_derivatives(a, false), arguments(O))
        O1 = operation(O)(args...)
        return simplify ? SymbolicUtils.simplify(O1) : O1
    end
end

function expand_derivatives(n::Num, simplify=false; occurrences=nothing)
    wrap(expand_derivatives(value(n), simplify; occurrences=occurrences))
end

function expand_derivatives(n::Complex{Num}, simplify=false; occurrences=nothing)
    wrap(ComplexTerm{Real}(expand_derivatives(real(n), simplify; occurrences=occurrences),
                           expand_derivatives(imag(n), simplify; occurrences=occurrences)))
end
_iszero(x) = false
_isone(x) = false

expand_derivatives(x, simplify=false;occurrences=nothing) = x

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
    istree(O) ? derivative(operation(O), (arguments(O)...,), Val(idx)) : 0
end

# Indicate that no derivative is defined.
struct NoDeriv
end
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
"""
function derivative(O, var; simplify=false)
    if O isa AbstractArray
        Num[Num(expand_derivatives(Differential(var)(value(o)), simplify)) for o in O]
    else
        Num(expand_derivatives(Differential(var)(value(O)), simplify))
    end
end

"""
$(SIGNATURES)

A helper function for computing the gradient of the expression `O` with respect to
an array of variable expressions.
"""
function gradient(O, vars::AbstractVector; simplify=false)
    Num[Num(expand_derivatives(Differential(v)(value(O)),simplify)) for v in vars]
end

"""
$(SIGNATURES)

A helper function for computing the Jacobian of an array of expressions with respect to
an array of variable expressions.
"""
function jacobian(ops::AbstractVector, vars::AbstractVector; simplify=false, scalarize=true)
    if scalarize
        ops = Symbolics.scalarize(ops)
        vars = Symbolics.scalarize(vars)
    end
    Num[Num(expand_derivatives(Differential(value(v))(value(O)),simplify)) for O in ops, v in vars]
end

function jacobian(ops, vars; simplify=false)
    ops = vec(scalarize(ops))
    vars = vec(scalarize(vars)) # Suboptimal, but prevents wrong results on Arr for now. Arr resulting from a symbolic function will fail on this due to unknown size.
    jacobian(ops, vars; simplify=simplify, scalarize=false)
end

"""
$(SIGNATURES)

A helper function for computing the sparse Jacobian of an array of expressions with respect to
an array of variable expressions.
"""
function sparsejacobian(ops::AbstractVector, vars::AbstractVector; simplify::Bool=false)
    ops = Symbolics.scalarize(ops)
    vars = Symbolics.scalarize(vars)
    sp = jacobian_sparsity(ops, vars)
    I,J,_ = findnz(sp)

    exprs = sparsejacobian_vals(ops, vars, I, J, simplify=simplify)

    sparse(I, J, exprs, length(ops), length(vars))
end

"""
$(SIGNATURES)

A helper function for computing the values of the sparse Jacobian of an array of expressions with respect to
an array of variable expressions given the sparsity structure.
"""
function sparsejacobian_vals(ops::AbstractVector, vars::AbstractVector, I::AbstractVector, J::AbstractVector; simplify::Bool=false)
    ops = Symbolics.scalarize(ops)
    vars = Symbolics.scalarize(vars)

    exprs = Num[]

    for (i,j) in zip(I, J)
        push!(exprs, Num(expand_derivatives(Differential(vars[j])(ops[i]), simplify)))
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


    simterm(x, f, args; kw...) = similarterm(x, f, args, symtype(x); kw...)

    # This rewriter notes down which u's appear in a
    # given du (whose index is stored in the `i` Ref)

    r = @rule ~x::(x->haskey(dict, x)) => begin
        push!(I, i[])
        push!(J, dict[~x])
        nothing
    end

    r =  Rewriters.Postwalk(r, similarterm=simterm)

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
"""
function hessian(O, vars::AbstractVector; simplify=false)
    vars = map(value, vars)
    first_derivs = map(value, vec(jacobian([values(O)], vars, simplify=simplify)))
    n = length(vars)
    H = Array{Num, 2}(undef,(n, n))
    fill!(H, 0)
    for i=1:n
        for j=1:i
            H[j, i] = H[i, j] = expand_derivatives(Differential(vars[i])(first_derivs[j]))
        end
    end
    H
end

isidx(x) = x isa TermCombination

basic_simterm(t, g, args; kws...) = Term{Any}(g, args)

let
    # we do this in a let block so that Revise works on the list of rules

    _scalar = one(TermCombination)

    linearity_rules = [
          @rule +(~~xs) => reduce(+, filter(isidx, ~~xs), init=_scalar)
          @rule *(~~xs) => reduce(*, filter(isidx, ~~xs), init=_scalar)
          @rule (~f)(~x::(!isidx)) => _scalar

          @rule (~f)(~x::isidx) => if haslinearity_1(~f)
              combine_terms_1(linearity_1(~f), ~x)
          else
              error("Function of unknown linearity used: ", ~f)
          end
          @rule (^)(~x::isidx, ~y) => ~y isa Number && isone(~y) ? ~x : (~x) * (~x)
          @rule (~f)(~x, ~y) => begin
              if haslinearity_2(~f)
                  a = isidx(~x) ? ~x : _scalar
                  b = isidx(~y) ? ~y : _scalar
                  combine_terms_2(linearity_2(~f), a, b)
              else
                  error("Function of unknown linearity used: ", ~f)
              end
          end
          @rule ~x::issym => 0]
    linearity_propagator = Fixpoint(Postwalk(Chain(linearity_rules); similarterm=basic_simterm))

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
        idx(i) = TermCombination(Set([Dict(i=>1)]))
        dict = Dict(u .=> idx.(1:length(u)))
        f = Rewriters.Prewalk(x->haskey(dict, x) ? dict[x] : x; similarterm=basic_simterm)(expr)
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
"""
function sparsehessian(op, vars::AbstractVector; simplify::Bool=false, full::Bool=true)
    op = value(op)
    vars = map(value, vars)
    S = hessian_sparsity(op, vars, full=full)
    I, J, _ = findnz(S)

    exprs = sparsehessian_vals(op, vars, I, J, simplify=simplify)

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
"""
function sparsehessian_vals(op, vars::AbstractVector, I::AbstractVector, J::AbstractVector; simplify::Bool=false)
    vars = Symbolics.scalarize(vars)

    exprs = Array{Num}(undef, length(I))
    fill!(exprs, 0)

    prev_j = 0
    d = nothing
    for (k, (i, j)) in enumerate(zip(I, J))
        j > i && continue
        if j != prev_j
            d = expand_derivatives(Differential(vars[j])(op), false)
        end
        expr = expand_derivatives(Differential(vars[i])(d), simplify)
        exprs[k] = expr
        prev_j = j
    end
    exprs
end

function SymbolicUtils.substitute(op::Differential, dict; kwargs...)
    @set! op.x = substitute(op.x, dict; kwargs...)
end
