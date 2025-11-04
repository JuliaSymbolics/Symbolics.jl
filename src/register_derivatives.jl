"""
    derivative_rule(::typeof(f), ::Val{NArgs}, args::ArgsT{VartypeT}, ::Val{I})

Define the derivative rule for `f` with `Nargs` arguments `args` with respect to the `I`th
argument. Do not define this function directly. Prefer using
[`@register_derivative`](@ref). Instead of calling this function directly, prefer
[`@derivative_rule`](@ref).
"""
function derivative_rule end

"""
    @register_derivative fn(args...) Ith_arg derivative

Register a symbolic derivative for a function. This typically accompanies a call to
[`@register_symbolic`](@ref) or [`@register_array_symbolic`](@ref) and defines how
[`expand_derivatives`](@ref) will behave when it tries to differentiate the registered
function.

The first argument to the macro is a call to the function whose derivative is being
defined. The call cannot have keyword arguments or default arguments. The call must have
either an exact number of arguments or a single variadic argument. For example, `f(a)`,
`f(a, b)`, `f(a, b, c)` and `f(args...)` are valid signatures. `f(a, b, args...)` is
invalid. If an exact number of arguments is provided, the defined derivative is specific
to that number of arguments. If the variadic signature is used, the defined derivative
is valid for all numbers of arguments. In case multiple derivatives are registered for
the same function, they must have different numbers of arguments. A derivative for an
exact number of arguments is more specific than a variadic definition. For example,
`@register_derivatives f(a, b) #...` is more specific than
`@register_derivatives f(args...) #...` for a 2-argument call to `f`. The arguments
can be referred to with their declared names inside the derivative definition.

The second argument to the macro is the argument with respect to which the derivative
rule is defined. For example, `@register_derivative f(a, b) 2 #...` is a derivative rule
with respect to the second argument of `f`. Mathematically, it represents
``\\frac{ \\partial f(a, b) }{ \\partial b }``. To define a generic derivative, this
argument can be an identifier. For example, `@register_derivative f(a, b) I #...` makes
`I` available in the derivative definition as the index of the argument with respect to
which the derivative is being taken.

The third argument to the macro is the derivative expression. This should be a symbolically
traceable expression returning the derivative of the specified function with respect to
the specified argument. In case of a variadic definition, the identifier `Nargs` is available
to denote the number of arguments provided to the function. In case the variadic form is
used, the arguments are available as a read-only array (mutation will error). Mutating
the array is unsafe and undefined behavior.

!!! note
    For functions that return arrays (such as those registered via `@register_array_symbolic`)
    the returned expression must be the Jacobian. Currently, support for differentiating array
    functions is considered experimental.

!!! warning
    The derivative expression MUST return a symbolic value, or `nothing` if the derivative is
    not defined. In case the result is a non-symbolic value, such as a constant derivative or
    Jacobian of array functions, the result MUST be wrapped in `Symbolics.SConst(..)`.

Following are example definitions of derivatives:

```julia
@register_derivative sin(x) 1 cos(x)
@register_derivative max(x, y) 2 ifelse(x >= y, 0, 1)
@register_derivative min(args...) I begin
  error("The rule for the derivative of `min` with \$Nargs arguments w.r.t the \$I-th argument is undefined.")
end
@register_derivative (foo::MyCallableStruct)(args...) I begin
  error("Oops! Didn't implement the derivative for \$foo")
end
```
"""
macro register_derivative(f::Expr, I::Union{Symbol, Int}, body)
    @assert Meta.isexpr(f, :call) """
    Incorrect `@register_derivative` syntax. The function must be provided as a call \
    signature. Got `$f` which is not a call signature.
    """
    fnhead = f.args[1]
    fncallargs = @view f.args[2:end]
    is_struct_der = Meta.isexpr(fnhead, :(::))
    if is_struct_der
        @assert length(fnhead.args) == 2 """
        Incorrect `@register_derivative` syntax. Registering derivatives of callable \
        structs requires providing a name for the struct. For example, instead of
        `@register_derivative (::MyStruct) # ...` use \
        `@register_derivative (x::MyStruct) # ...`.
        """
    end
    @assert !any(Base.Fix2(Meta.isexpr, :kw), fncallargs) """
    Incorrect `@register_derivative` syntax. The function cannot have default arguments.
    """
    @assert !Meta.isexpr(fncallargs[1], :parameters) """
    Incorrect `@register_derivative` syntax. The function cannot have keyword arguments.
    """

    is_varargs = Meta.isexpr(fncallargs[end], :...)
    if is_varargs
        @assert length(fncallargs) == 1 """
        Incorrect `@register_derivative` syntax. The function call signature must either \
        be a single variadic argument `@register_derivative foo(args...) #...` or a \
        concrete number of arguments `@register_derivative foo(arg1, arg2, arg3) # ...`.
        """
    end

    derhead = Expr(:call, :($Symbolics.derivative_rule), is_struct_der ? fnhead : :(::($typeof($fnhead))))
    Nargs = is_varargs ? :Nargs : length(fncallargs)
    push!(derhead.args, :(::Val{$Nargs}))
    args_name = gensym(:args)
    push!(derhead.args, :($args_name::$SymbolicUtils.ROArgsT{$VartypeT}))
    push!(derhead.args, :(::Val{$I}))

    if is_varargs || I isa Symbol
        derhead = Expr(:where, derhead)
        is_varargs && push!(derhead.args, Nargs)
        I isa Symbol && push!(derhead.args, I)
    end

    if is_varargs
        unpack = Expr(:(=), fncallargs[1].args[1], args_name)
    else
        unpack = Expr(:tuple)
        append!(unpack.args, fncallargs)
        unpack = Expr(:(=), unpack, args_name)
    end

    return esc(Expr(:function, derhead, Expr(:block, unpack, body)))
end

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

        # Using the macro here doesn't work somehow.
        @eval function derivative_rule(::typeof($modu.$fun), ::Val{$arity}, args::SymbolicUtils.ArgsT{VartypeT}, ::Val{$i})
            $SConst($expr)
        end
    end
end

Base.@propagate_inbounds function _derivative_rule_proxy(f, args::NTuple{N, SymbolicT}, ::Val{I}) where {N, I}
    _derivative_rule_proxy(f, Val{N}(), args, Val{I}())
end
Base.@propagate_inbounds function _derivative_rule_proxy(f, ::Val{N}, args::NTuple{N, SymbolicT}, ::Val{I}) where {N, I}
    _derivative_rule_proxy(f, Val{N}(), SymbolicUtils.ArgsT{VartypeT}(args), Val{I}())
end
Base.@propagate_inbounds function _derivative_rule_proxy(f, args::Tuple, ::Val{I}) where {I}
    _derivative_rule_proxy(f, Val{length(args)}(), args, Val{I}())
end
Base.@propagate_inbounds function _derivative_rule_proxy(f, ::Val{N}, args::Tuple{Vararg{Any, N}}, ::Val{I}) where {N, I}
    args = ntuple(BSImpl.Const{VartypeT} ∘ Base.Fix1(getindex, args), Val{N}())
    _derivative_rule_proxy(f, Val{N}(), args, Val{I}())
end
Base.@propagate_inbounds function _derivative_rule_proxy(f, args::ROArgsT{VartypeT}, ::Val{I}) where {I}
    @inbounds _derivative_rule_proxy(f, Val{length(args)}(), args, Val{I}())
end
Base.@propagate_inbounds function _derivative_rule_proxy(f, ::Val{N}, args::ROArgsT{VartypeT}, ::Val{I}) where {N, I}
    @boundscheck checkbounds(args, N)
    derivative_rule(f, Val{N}(), args, Val{I}())
end
Base.@propagate_inbounds function _derivative_rule_proxy(f, args::ArgsT{VartypeT}, ::Val{I}) where {I}
    @inbounds _derivative_rule_proxy(f, Val{length(args)}(), args, Val{I}())
end
Base.@propagate_inbounds function _derivative_rule_proxy(f, ::Val{N}, args::ArgsT{VartypeT}, ::Val{I}) where {N, I}
    @boundscheck checkbounds(args, N)
    _derivative_rule_proxy(f, Val{N}(), ROArgsT{VartypeT}(args), Val{I}())
end
Base.@propagate_inbounds function _derivative_rule_proxy(f, args::AbstractArray{SymbolicT}, ::Val{I}) where {I}
    @inbounds _derivative_rule_proxy(f, Val{length(args)}(), args, Val{I}())
end
Base.@propagate_inbounds function _derivative_rule_proxy(f, ::Val{N}, args::AbstractArray{SymbolicT}, ::Val{I}) where {N, I}
    @boundscheck checkbounds(args, N)
    _derivative_rule_proxy(f, Val{N}(), ArgsT{VartypeT}(args), Val{I}())
end
Base.@propagate_inbounds function _derivative_rule_proxy(f, args::AbstractArray, ::Val{I}) where {I}
    @inbounds _derivative_rule_proxy(f, Val{length(args)}(), args, Val{I}())
end
Base.@propagate_inbounds function _derivative_rule_proxy(f, ::Val{N}, args::AbstractArray, ::Val{I}) where {N, I}
    @boundscheck checkbounds(args, N)
    _args = ArgsT{VartypeT}()
    sizehint!(_args, N)
    for a in args
        push!(_args, BSImpl.Const{VartypeT}(a))
    end
    _derivative_rule_proxy(f, Val{N}(), _args, Val{I}())
end

"""
    @derivative_rule f(args...) I

Query Symbolics.jl's derivative rule system for the derivative of `f(args...)` with respect to
`args[I]`. Returns a symbolic result representing the derivative. In case the derivative rule is
not defined, evaluates to `nothing`.

The first argument to the macro must be a valid function call syntax. Splatting of arguments is
permitted. The second argument must be an expression or literal evaluating to the index of the
argument with respect to which the derivative is required.

The derivative rule can dispatch statically if `f`, the number of arguments and `I` are known
at compile time. Example invocations are:

```julia
# static dispatch
@derivative_rule sin(x) 1
# static dispatch if `xs` is a tuple
@derivative_rule max(xs...) 2
# static dispatch if `y` and `w` are tuples, and `N + 2K` is a compile-time constant
@derivative_rule foo(x, y..., z, w...) (N + 2K)
```
"""
macro derivative_rule(f, I)
    @assert Meta.isexpr(f, :call) """
    Incorrect `@derivative_rule` syntax. The function must be provided as a call \
    signature. Got `$f` which is not a call signature.
    """
    fnhead = f.args[1]
    fncallargs = @view f.args[2:end]
    result = Expr(:call, _derivative_rule_proxy, fnhead)
    if length(fncallargs) == 1 && Meta.isexpr(fncallargs[1], :...)
        push!(result.args, fncallargs[1].args[1])
    elseif any(Base.Fix2(Meta.isexpr, :...), fncallargs)
        args = Expr(:tuple)
        append!(args.args, fncallargs)
        push!(result.args, args)
    else
        push!(result.args, :(Val{$(length(fncallargs))}()))
        args = Expr(:tuple)
        append!(args.args, fncallargs)
        push!(result.args, args)
    end
    push!(result.args, :(Val{$I}()))
    return esc(result)
end

@register_derivative +(args...) I COMMON_ONE
@register_derivative *(args...) I begin
    if I == 1
        SymbolicUtils.mul_worker(VartypeT, view(args, 2:Nargs))
    elseif I == Nargs
        SymbolicUtils.mul_worker(VartypeT, view(args, 1:(Nargs-1)))
    else
        t1 = SymbolicUtils.mul_worker(VartypeT, view(args, 1:(Nargs-1)))
        t2 = SymbolicUtils.mul_worker(VartypeT, view(args, 2:Nargs))
        t1 * t2
    end
end
@register_derivative one(x) 1 COMMON_ZERO

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
@inline derivative_idx(::Any, ::Any) = COMMON_ZERO
function derivative_idx(O::VartypeT, idx::Int)
    iscall(O) || return COMMON_ZERO
    f = operation(O)
    args = arguments(O)
    return @derivative_rule f(args...) idx
end

