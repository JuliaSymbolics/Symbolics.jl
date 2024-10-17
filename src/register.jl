using SymbolicUtils: Symbolic

"""
    @register_symbolic(expr, define_promotion = true, Ts = [Real])

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
See `@register_array_symbolic` to register functions which return arrays.
"""
macro register_symbolic(expr, define_promotion = true, Ts = :([]), wrap_arrays = true)
    f, ftype, argnames, Ts, ret_type = destructure_registration_expr(expr, Ts)

    args′ = map((a, T) -> :($a::$T), argnames, Ts)
    ret_type = isnothing(ret_type) ? Real : ret_type

    fexpr = :(Symbolics.@wrapped function $f($(args′...))
                  args = [$(argnames...),]
                  unwrapped_args = map($nested_unwrap, args)
                  res = if !any($is_symbolic_or_array_of_symbolic, unwrapped_args)
                      $f(unwrapped_args...) # partial-eval if all args are unwrapped
                  else
                      $Term{$ret_type}($f, unwrapped_args)
                  end
                  if typeof.(args) == typeof.(unwrapped_args)
                      return res
                  else
                      return $wrap(res)
                  end
              end $wrap_arrays)

    if define_promotion
        fexpr = :($fexpr; (::$typeof($promote_symtype))(::$ftype, args...) = $ret_type)
    end
    esc(fexpr)
end

function destructure_registration_expr(expr, Ts)
    if expr.head === :(::)
        ret_type = expr.args[2]
        expr = expr.args[1]
    else
        ret_type = nothing
    end
    @assert expr.head === :call
    @assert Ts.head === :vect
    Ts = Ts.args

    f = expr.args[1]
    args = expr.args[2:end]

    # Default arg types to Real
    Ts = map(a -> a isa Symbol ? Real : (@assert(a.head == :(::)); a.args[2]), args)
    argnames = map(a -> a isa Symbol ? a : a.args[1], args)

    ftype = if f isa Expr && f.head == :(::)
        if length(f.args) == 1
            error("please name the callable object, i.e. use (f::$(f.args[end])) instead of $f")
        end
        @assert length(f.args) == 2
        f.args[end]
    else
        :($typeof($f))
    end
    f, ftype, argnames, Ts, ret_type
end

nested_unwrap(x) = unwrap(x)
nested_unwrap(x::Arr) = unwrap(x)
nested_unwrap(x::AbstractArray) = unwrap.(x)

function is_symbolic_or_array_of_symbolic(x)
    return issym(x) || iscall(x)
end
function is_symbolic_or_array_of_symbolic(arr::AbstractArray)
    return any(is_symbolic_or_array_of_symbolic.(arr))
end

symbolic_eltype(x) = eltype(x)
symbolic_eltype(::AbstractArray{symT}) where {eT, symT <: SymbolicUtils.Symbolic{eT}} = eT
symbolic_eltype(::AbstractArray{Num}) = Real
symbolic_eltype(::AbstractArray{symT}) where {eT, symT <: Arr{eT}} = eT

function register_array_symbolic(f, ftype, argnames, Ts, ret_type, partial_defs = :(), define_promotion = true, wrap_arrays = true)
    def_assignments = MacroTools.rmlines(partial_defs).args
    defs = map(def_assignments) do ex
        @assert ex.head == :(=)
        ex.args[1] => ex.args[2]
    end |> Dict


    args′ = map((a, T) -> :($a::$T), argnames, Ts)
    fexpr = quote
        @wrapped function $f($(args′...))
            args = [$(argnames...),]
            unwrapped_args = map($nested_unwrap, args)
            eltype = $symbolic_eltype
            res = if !any($is_symbolic_or_array_of_symbolic, unwrapped_args)
                $f(unwrapped_args...) # partial-eval if all args are unwrapped
            elseif $ret_type == nothing || ($ret_type <: AbstractArray)
                $array_term($(Expr(:parameters, [Expr(:kw, k, v) for (k, v) in defs]...)), $f, unwrapped_args...)
            else
                $Term{$ret_type}($f, unwrapped_args)
            end

            if typeof.(args) == typeof.(unwrapped_args)
                return res
            else
                return $wrap(res)
            end
        end $wrap_arrays
    end |> esc

    if define_promotion
        container_type = get(defs, :container_type, :($propagate_atype(f, $(argnames...))))
        etype = get(defs, :eltype, :($propagate_eltype(f, $(argnames...))))
        ndim = get(defs, :ndims, nothing)
        is_callable_struct = f isa Expr && f.head == :(::)
        fn_arg = if is_callable_struct
            f
        else
            :(f::$ftype)
        end
        fn_arg_name = if is_callable_struct
            f.args[1]
        else
            :f
        end
        promote_expr = quote
            function (::$typeof($promote_symtype))($fn_arg, $(argnames...))
                f = $fn_arg_name
                container_type = $container_type
                etype = $etype
                $(
                    if ndim === nothing
                        :(return container_type{etype})
                    else
                        :(ndim = $ndim; return container_type{etype, ndim})
                    end
                )
            end
        end |> esc
        fexpr = :($fexpr; $promote_expr)
    end

    return fexpr
end

"""
    @register_array_symbolic(expr, define_promotion = true)

Example:

```julia
# Let's say vandermonde takes an n-vector and returns an n x n matrix
@register_array_symbolic vandermonde(x::AbstractVector) begin
    size=(length(x), length(x))
    eltype=eltype(x) # optional, will default to the promoted eltypes of x
end
```

You can also register calls on callable structs:

```julia
@register_array_symbolic (c::Conv)(x::AbstractMatrix) begin
    size=size(x) .- size(c.kernel) .+ 1
    eltype=promote_type(eltype(x), eltype(c))
end
```

If `define_promotion = true` then a promotion method in the form of
```julia
SymbolicUtils.promote_symtype(::typeof(f_registered), args...) = # inferred or annotated return type
```

is defined for the register function. Note that when defining multiple register
overloads for one function, all the rest of the registers must set
`define_promotion` to `false` except for the first one, to avoid method
overwriting.
"""
macro register_array_symbolic(expr, block, define_promotion = true, wrap_arrays = true)
    f, ftype, argnames, Ts, ret_type = destructure_registration_expr(expr, :([]))
    register_array_symbolic(f, ftype, argnames, Ts, ret_type, block, define_promotion, wrap_arrays)
end
