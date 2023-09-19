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
macro register_symbolic(expr, defs = true, flags...)
    if expr.head === :(::)
        ret_type = expr.args[2]
        expr = expr.args[1]
    else
        ret_type = nothing
    end
    @assert expr.head === :call

    f = expr.args[1]
    args = expr.args[2:end]

    # Default arg types to Real
    Ts = map(a -> a isa Symbol ? Real : (@assert(a.head == :(::)); a.args[2]), args)
    argnames = map(a -> a isa Symbol ? a : a.args[1], args)


    ftype = if f isa Expr && f.head == :(::)
        @assert length(f.args) == 2
        f.args[end]
    else
        :($typeof($f))
    end

    if :array in flags
        return register_array_symbolic(f, ftype, argnames, Ts, ret_type, defs)
    end

    args′ = map((a, T) -> :($a::$T), argnames, Ts)
    ret_type = isnothing(ret_type) ? Real : ret_type

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

    if defs
        fexpr = :($fexpr; (::$typeof($promote_symtype))(::$ftype, args...) = $ret_type)
    end
    esc(fexpr)
end

function register_array_symbolic(f, ftype, argnames, Ts, ret_type, partial_defs = :())
    def_assignments = MacroTools.rmlines(partial_defs).args
    defs = map(def_assignments) do ex
        @assert ex.head == :(=)
        ex.args[1] => ex.args[2]
    end |> Dict

    if haskey(defs, :size)
        # we don't store size, instead we store extents of indices "shape"
        # or axes -- but it can also be Unknown()
        shape = quote
            Tuple(map(x->1:x, Symbolics.@oops $(defs[:size])))
        end
    elseif haskey(defs, :shape)
        shape = defs[:shape]
    else
        shape = Unknown()
    end

    eltype = get(defs, :eltype, Unknown())

    ndims = get(defs, :ndims, Unknown())

    atype = get(defs, :container_type, Unknown())

    args′ = map((a, T) -> :($a::$T), argnames, Ts)
    quote
        @wrapped function $f($(args′...))
            args = [$(argnames...),]
            unwrapped_args = map($unwrap, args)
            res = if !any(x->$issym(x) || $istree(x), unwrapped_args)
                $f(unwrapped_args...) # partial-eval if all args are unwrapped
            elseif $ret_type == nothing || ($ret_type <: AbstractArray)
                $arrterm($f, unwrapped_args...; atype=$atype, eltype=$eltype, ndims=$ndims, shape=$shape)
            else
                $Term{$ret_type}($f, unwrapped_args)
            end

            if typeof.(args) == typeof.(unwrapped_args)
                return res
            else
                return $wrap(res)
            end
        end
    end |> esc
end

Base.@deprecate_binding var"@register" var"@register_symbolic"
