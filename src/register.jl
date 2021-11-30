using SymbolicUtils: Symbolic

"""
    @register(expr, define_promotion = true, Ts = [Num, Symbolic, Real])

Overload appropriate methods so that Symbolics can stop tracing into the
registered function. If `define_promotion` is true, then a promotion method in
the form of
```julia
SymbolicUtils.promote_symtype(::typeof(f_registered), args...) = Real # or the annotated return type
```
is defined for the register function. Note that when defining multiple register
overloads for one function, all the rest of the registers must set
`define_promotion` to `false` except for the first one, to avoid method
overwritting.

# Examples
```julia
@register foo(x, y)
@register foo(x, y::Bool) false # do not overload a duplicate promotion rule
@register goo(x, y::Int) # `y` is not overloaded to take symbolic objects
@register hoo(x, y)::Int # `hoo` returns `Int`
```
"""
macro register(expr, define_promotion = true, Ts = [])
    if expr.head === :(::)
        ret_type = expr.args[2]
        expr = expr.args[1]
    else
        ret_type = Real
    end

    @assert expr.head === :call

    f = expr.args[1]
    args = expr.args[2:end]

    types = map(args) do x
        if x isa Symbol
            :(($Real, $wrapper_type($Real), $Symbolic{<:$Real}))
        elseif Meta.isexpr(x, :(::))
            T = x.args[2]
            :($has_symwrapper($T) ?
              ($T, $Symbolic{<:$T}, $wrapper_type($T),) :
              ($T, $Symbolic{<:$T}))
        else
            error("Invalid argument format $x")
        end
    end

    eval_method = :(@eval function $f($(Expr(:$, :(s...))),)
                        $wrap($Term{$ret_type}($f, [$(Expr(:$, :(s_syms...)))]))
                    end)
    verbose = false
    mod, fname = f isa Expr && f.head == :(.) ? f.args : (:(@__MODULE__), QuoteNode(f))
    quote
        __Ts = [Tuple{x...} for x in Iterators.product($(types...),)
                if !all(x->Symbolics.has_symwrapper(x) &&
                        Symbolics.wrapper_type(x) != x, x)]
        if $verbose
            println("Candidates")
            map(println, __Ts)
        end

        for sig in __Ts
            s = map(((i,T,),)->Expr(:(::), Symbol("arg", i), T), enumerate(sig.parameters))
            s_syms = map(x->x.args[1], s)
            $eval_method
        end
        if $define_promotion
            (::$typeof($promote_symtype))(::$typeof($f), args...) = $ret_type
        end
    end |> esc
end

# Ensure that Num that get @registered from outside the ModelingToolkit
# module can work without having to bring in the associated function into the
# ModelingToolkit namespace. We basically store information about functions
# registered at runtime in a ModelingToolkit variable,
# `registered_external_functions`. It's not pretty, but we are limited by the
# way GeneralizedGenerated builds a function (adding "ModelingToolkit" to every
# function call).
# ---
const registered_external_functions = Dict{Symbol,Module}()
function inject_registered_module_functions(expr)
    MacroTools.postwalk(expr) do x
        # Find all function calls in the expression and extract the function
        # name and calling module.
        MacroTools.@capture(x, f_module_.f_name_(xs__))
        if isnothing(f_module)
            MacroTools.@capture(x, f_name_(xs__))
        end

        if !isnothing(f_name)
            # Set the calling module to the module that registered it.
            mod = get(registered_external_functions, f_name, f_module)
            if !isnothing(mod) && mod != Base
                x.args[1] = :($mod.$f_name)
            end
        end

        return x
    end
end
