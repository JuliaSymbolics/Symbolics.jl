using SymbolicUtils: Symbolic

"""
    @register_symbolic(expr, define_promotion = true)

Overload appropriate methods so that Symbolics can stop tracing into the
registered function. If `define_promotion` is true, then a promotion method in
the form of
```julia
SymbolicUtils.promote_symtype(::typeof(f_registered), args...) = Real # or the annotated return type
```
is defined for the register function. Note that when defining multiple register
overloads for one function, all the rest of the registers must set
`define_promotion` to `false` except for the first one, to avoid method
overwriting.

# Examples
```julia
@register_symbolic foo(x, y)
@register_symbolic foo(x, y::Bool) false # do not overload a duplicate promotion rule
@register_symbolic goo(x, y::Int) # `y` is not overloaded to take symbolic objects
@register_symbolic hoo(x, y)::Int # `hoo` returns `Int`
```
"""
macro register_symbolic(expr, define_promotion = true)
    if expr.head === :(::)
        ret_type = expr.args[2]
        expr = expr.args[1]
    else
        ret_type = Real
    end
    @assert expr.head === :call

    f = expr.args[1]
    args = expr.args[2:end]

    # Default arg types to Real
    Ts = map(a -> a isa Symbol ? Real : (@assert(a.head == :(::)); a.args[2]), args)
    argnames = map(a -> a isa Symbol ? a : a.args[1], args)
    args′ = map((a, T) -> :($a::$T), argnames, Ts)

    # @register_symbolic f(x::T1, y::T2)
    #
    # `types`: for every argument find the types that
    # could be taken as arguments. These are:
    #
    # 1) T 2) wrapper_type(T) 3) Symbolic{T}
    #
    # However later while emiting methods we omit the ones
    # that are all 1) since those are expected to be defined
    # outside Symbolics

    ftype = if f isa Expr && f.head == :(::)
        @assert length(f.args) == 2
        f.args[end]
    else
        :($typeof($f))
    end
    fexpr = :(@wrapped $f($(args′...)) = $Term{$ret_type}($f, [$(argnames...)]))

    if define_promotion
        fexpr = :($fexpr; (::$typeof($promote_symtype))(::$ftype, args...) = $ret_type)
    end
    esc(fexpr)
end

Base.@deprecate_binding var"@register" var"@register_symbolic"

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
