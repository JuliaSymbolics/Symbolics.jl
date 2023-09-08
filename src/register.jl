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


    ftype = if f isa Expr && f.head == :(::)
        @assert length(f.args) == 2
        f.args[end]
    else
        :($typeof($f))
    end
    fexpr = :(@wrapped function $f($(args′...))
                  args = [$(argnames...),]
                  unwrapped_args = map($unwrap, args)
                  res = if !any(x->$issym(x) || $istree(x), unwrapped_args)
                      $f(unwrapped_args...) # partial-eval if all args are unwrapped
                  else
                      $Term{$ret_type}($f, unwrapped_args)
                  end
                  if typeof.(args) == typeof.(unwrapped_args)
                      return res
                  else
                      return $wrap(res)
                  end
              end)

    if define_promotion
        fexpr = :($fexpr; (::$typeof($promote_symtype))(::$ftype, args...) = $ret_type)
    end
    esc(fexpr)
end

Base.@deprecate_binding var"@register" var"@register_symbolic"

