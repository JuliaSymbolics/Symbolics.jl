export @symbolic_wrap, @wrapped

# Turn A{X} <: B{Int, X} into
#
# B{Int, X} where X
function set_where(subt, supert)
    if !(supert isa Expr && supert.head == :curly)
        return supert
    end

    Ss = []
    if subt isa Expr && subt.head == :curly
        Ss = subt.args[2:end]
    end

    Ts = intersect(supert.args[2:end], Ss)
    Expr(:where, supert, Ts...)
end

getname(x::Symbol) = x

function getname(x::Expr)
    @assert x.head == :curly
    return x.args[1]
end

macro symbolic_wrap(expr)
    @assert expr isa Expr && expr.head == :struct
    @assert expr.args[2].head == :(<:)
    sig = expr.args[2]
    T = getname(sig.args[1])
    supertype = set_where(sig.args[1], sig.args[2])

    quote
        $expr

        Symbolics.has_symwrapper(::Type{<:$supertype}) = true
        Symbolics.wrapper_type(::Type{<:$supertype}) = $T
        Symbolics.is_wrapper_type(::Type{<:$T}) = true # used in `@register`
        Symbolics.wraps_type(::Type{$T}) = $supertype
        Symbolics.iswrapped(::$T) = true
    end |> esc
end

iswrapped(x) = false
unwrap(x) = x

function wrap(x)
    T = SymbolicUtils.symtype(x)
    Symbolics.has_symwrapper(T) ?
        Symbolics.wrapper_type(T)(x) : x
end

function wrapper_type end
function wraps_type end

has_symwrapper(::Type) = false
is_wrapper_type(::Type) = false

# helper function to extract keyword argument names from expressions
function expr_kwargname(kwarg)
    if kwarg isa Expr && kwarg.head == :kw
        kwarg.args[1]
    elseif kwarg isa Expr && kwarg.head == :(...)
        kwarg.args[1]
    else
        kwarg
    end
end

# helper function to extract argument names from expressions
function expr_argname(arg)
    if arg isa Expr && (arg.head == :(::) || arg.head == :(...))
        arg.args[1]
    elseif arg isa Expr
        error("$arg not supported as an argument")
    else
        arg
    end
end

function wrap_func_expr(mod, expr)
    @assert expr.head == :function || (expr.head == :(=) &&
                                       expr.args[1] isa Expr &&
                                       expr.args[1].head == :call)

    def = splitdef(expr)
    body = def[:body]
    fname = def[:name]
    args = get(def, :args, [])
    kwargs = get(def, :kwargs, [])
    args_names = expr_argname.(args)
    kwargs_names = expr_kwargname.(kwargs)
    
    wrap_func_expr(mod, fname, args, kwargs, args_names, kwargs_names, Symbol[], body)
end

function wrap_func_expr(
    mod, fname, args, kwargs, args_names, kwargs_names, whereparams, body;
    abstract_arg_types=nothing
)
    names = vcat(args_names, kwargs_names)
    
    function type_options(wparams, arg, arg_ind)
        pmod = parentmodule(mod)
        atype = isnothing(abstract_arg_types) ? Any : abstract_arg_types[arg_ind]
        if arg isa Expr && arg.head == :(::)
            T = Base.eval(mod, quote
                let $(Symbol(pmod)) = $(pmod); # make name of parent module available in eval scope   
                    #=
                    NOTE
                    `typeintersect` is important here for consecutive calls to `specialize_methods` 
                    with conceptually different super types. 
                    E.g.: Consider we first specialize `*(::AbstractMatrix, ::AbstractVector)` to 
                    redirect to `_matvec`, and then `*(::AbstractMatrix, ::AbstractMatrix)` to 
                    redirect to `_matmul`. If we encounter some existing method for `*` which accepts
                    an `AbstractMatrix` and an `VecOrMat` (type union), then we accidentally redirect 
                    a matrix-vector-product to `_matmul` without `typeintersect`.
                    =#
                    typeintersect($(atype), $(arg.args[2]) where {$(wparams...)})
                end
            end)
            has_symwrapper(T) ? (T, :(SymbolicUtils.Symbolic{<:$T}), wrapper_type(T)) :
                                (T, :(SymbolicUtils.Symbolic{<:$T}))
        elseif arg isa Expr && arg.head == :(...)
            Ts = type_options(wparams, arg.args[1], arg_ind)
            map(x->Vararg{x},Ts)
        else
            (Any,)
        end
    end

    types = [type_options(whereparams, arg, arg_ind) for (arg_ind, arg)=enumerate(args)]

    impl_name = Symbol(fname,"_", hash(string(args)*string(kwargs)*string(types)))

    impl = :(function $impl_name($(names...))
        $body
    end)
    # TODO: maybe don't drop first lol
    methods = map(Iterators.drop(Iterators.product(types...), 1)) do Ts
        method_args = map(names, Ts) do n, T
            :($n::$T)
        end

        fbody = :(if any($iswrapped, ($(names...),))
                      $wrap($impl_name($([:($unwrap($arg)) for arg in names]...)))
                  else
                      $impl_name($(names...))
                  end)

        if isempty(kwargs)
            :(function $fname($(method_args...))
                  $fbody
              end)
        else
            :(function $fname($(method_args...); $(kwargs...))
                  $fbody
              end)
        end
    end

    quote
        $impl
        $(methods...)
    end
end

macro wrapped(expr)
    esc(wrap_func_expr(__module__, expr))
end
