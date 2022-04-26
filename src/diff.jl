abstract type Operator <: Function end

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
(D::Differential)(x) = Term{symtype(x)}(D, [x])
(D::Differential)(x::Num) = Num(D(value(x)))
SymbolicUtils.promote_symtype(::Differential, x) = x

is_derivative(x::Term) = operation(x) isa Differential
is_derivative(x) = false

Base.:*(D1, D2::Differential) = D1 ∘ D2
Base.:*(D1::Differential, D2) = D1 ∘ D2
Base.:*(D1::Differential, D2::Differential) = D1 ∘ D2
Base.:^(D::Differential, n::Integer) = _repeat_apply(D, n)

Base.show(io::IO, D::Differential) = print(io, "Differential(", D.x, ")")

Base.:(==)(D1::Differential, D2::Differential) = isequal(D1.x, D2.x)
Base.hash(D::Differential, u::UInt) = hash(D.x, xor(u, 0xdddddddddddddddd))

_isfalse(occ::Bool) = occ === false
_isfalse(occ::Term) = _isfalse(operation(occ))

function occursin_info(x, expr)
    if symtype(expr) <: AbstractArray
        error("Differentiation of expressions involving arrays and array variables is not yet supported.")
    end

    # Allow scalarized expressions
    function is_scalar_indexed(ex)
        (istree(ex) && operation(ex) == getindex && !(symtype(ex) <: AbstractArray)) ||
        (istree(ex) && (issym(operation(ex)) || istree(operation(ex))) &&
         is_scalar_indexed(operation(ex)))
    end

    if is_scalar_indexed(x) && is_scalar_indexed(expr) &&
        isequal(first(arguments(x)), first(arguments(expr)))
        return isequal(arguments(x), arguments(expr))
    end
    if is_scalar_indexed(x) && is_scalar_indexed(expr) &&
        !occursin(first(arguments(x)), first(arguments(expr)))
        return false
    end

    if is_scalar_indexed(expr) && !is_scalar_indexed(x) && !occursin(x, expr)
        return false
    end

    !istree(expr) && return false
    if isequal(x, expr)
        true
    else
        args = map(a->occursin_info(x, a), arguments(expr))
        if all(_isfalse, args)
            return false
        end
        Term{Real}(true, args)
    end
end

function occursin_info(x, expr::Sym)
    if symtype(expr) <: AbstractArray
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
        if O isa Union{Add, Mul}
            any(recursive_hasoperator(op), keys(O.dict))
        elseif O isa Pow
            recursive_hasoperator(op)(O.base) || recursive_hasoperator(op)(O.exp)
        elseif O isa SymbolicUtils.Div
            recursive_hasoperator(op)(O.num) || recursive_hasoperator(op)(O.den)
        else
            any(recursive_hasoperator(op), arguments(O))
        end
    end
end

"""
$(SIGNATURES)

TODO
"""
function expand_derivatives(O::Symbolic, simplify=false; occurances=nothing)
    if istree(O) && isa(operation(O), Differential)
        arg = only(arguments(O))
        arg = expand_derivatives(arg, false)

        if occurances == nothing
            occurances = occursin_info(operation(O).x, arg)
        end

        _isfalse(occurances) && return 0
        occurances isa Bool && return 1 # means it's a `true`

        D = operation(O)

        if !istree(arg)
            return D(arg) # Cannot expand
        elseif (op = operation(arg); isa(op, Sym))
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
            return expand_derivatives(O, simplify; occurances)
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
            t2 = expand_derivatives(D(inner_args[i]),false, occurances=arguments(occurances)[i])

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

function expand_derivatives(n::Num, simplify=false; occurances=nothing)
    Num(expand_derivatives(value(n), simplify; occurances=occurances))
end

_iszero(x) = false
_isone(x) = false

expand_derivatives(x, simplify=false;occurances=nothing) = x

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
for (modu, fun, arity) ∈ DiffRules.diffrules()
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

_repeat_apply(f, n) = n == 1 ? f : f ∘ _repeat_apply(f, n-1)
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

A helper function for computing the derivative of an expression with respect to
`var`.
"""
function derivative(O, v; simplify=false)
    if O isa AbstractArray
        Num[Num(expand_derivatives(Differential(v)(value(o)), simplify)) for o in O]
    else
        Num(expand_derivatives(Differential(v)(value(O)), simplify))
    end
end

"""
$(SIGNATURES)

A helper function for computing the gradient of an expression with respect to
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
function jacobian(ops::AbstractVector, vars::AbstractVector; simplify=false)
    ops = Symbolics.scalarize(ops)
    vars = Symbolics.scalarize(vars)
    Num[Num(expand_derivatives(Differential(value(v))(value(O)),simplify)) for O in ops, v in vars]
end

function jacobian(ops::ArrayLike{T, 1}, vars::ArrayLike{T, 1}; simplify=false) where T
    ops = scalarize(ops)
    vars = scalarize(vars) # Suboptimal, but prevents wrong results on Arr for now. Arr resulting from a symbolic function will fail on this due to unknown size.
    Num[Num(expand_derivatives(Differential(value(v))(value(O)),simplify)) for O in ops, v in vars]
end

"""
$(SIGNATURES)

A helper function for computing the sparse Jacobian of an array of expressions with respect to
an array of variable expressions.
"""
function sparsejacobian(ops::AbstractVector, vars::AbstractVector; simplify=false)
    I = Int[]
    J = Int[]
    du = Num[]

    ops = Symbolics.scalarize(ops)
    vars = Symbolics.scalarize(vars)
    sp = jacobian_sparsity(ops, vars)
    I,J,_ = findnz(sp)

    exprs = Num[]

    for (i,j) in zip(I, J)
        push!(exprs, Num(expand_derivatives(Differential(vars[j])(ops[i]), simplify)))
    end
    sparse(I, J, exprs, length(ops), length(vars))
end

"""
```julia
jacobian_sparsity(ops::AbstractVector, vars::AbstractVector)
```

Return the sparsity pattern of the Jacobian of an array of expressions with respect to
an array of variable expressions.
"""
function jacobian_sparsity(du, u)
    du = map(value, du)
    u = map(value, u)
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
```julia
jacobian_sparsity(op!,output::Array{T},input::Array{T}) where T<:Number
```

Return the sparsity pattern of the Jacobian of the mutating function `op!(output,input,args...)`.
"""
function jacobian_sparsity(op!,output::Array{T},input::Array{T}, args...) where T<:Number
    eqs=similar(output,Num)
    fill!(eqs,false)
    vars=ArrayInterface.restructure(input,[variable(i) for i in eachindex(input)])
    op!(eqs,vars, args...)
    jacobian_sparsity(eqs,vars)
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

A helper function for computing the Hessian of an expression with respect to
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

"""
    hessian_sparsity(ops::AbstractVector, vars::AbstractVector)

Return the sparsity pattern of the Hessian of an array of expressions with respect to
an array of variable expressions.
"""
function hessian_sparsity end
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
          @rule ~x::(x->x isa Sym) => 0]
    linearity_propagator = Fixpoint(Postwalk(Chain(linearity_rules); similarterm=basic_simterm))

    global hessian_sparsity

    function hessian_sparsity(f, u)
        @assert !(f isa AbstractArray)
        f = value(f)
        u = map(value, u)
        idx(i) = TermCombination(Set([Dict(i=>1)]))
        dict = Dict(u .=> idx.(1:length(u)))
        f = Rewriters.Prewalk(x->haskey(dict, x) ? dict[x] : x; similarterm=basic_simterm)(f)
        lp = linearity_propagator(f)
        _sparse(lp, length(u))
    end
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
function sparsehessian(O, vars::AbstractVector; simplify=false)
    O = value(O)
    vars = map(value, vars)
    S = hessian_sparsity(O, vars)
    I, J, _ = findnz(S)
    exprs = Array{Num}(undef, length(I))
    fill!(exprs, 0)
    prev_j = 0
    d = nothing
    for (k, (i, j)) in enumerate(zip(I, J))
        j > i && continue
        if j != prev_j
            d = expand_derivatives(Differential(vars[j])(O), false)
        end
        expr = expand_derivatives(Differential(vars[i])(d), simplify)
        exprs[k] = expr
        prev_j = j
    end
    H = sparse(I, J, exprs, length(vars), length(vars))
    for (i, j) in zip(I, J)
        j > i && (H[i, j] = H[j, i])
    end
    return H
end
