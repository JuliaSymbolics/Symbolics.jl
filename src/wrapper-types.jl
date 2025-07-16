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

function SymbolicIndexingInterface.getname(x::Expr)
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

function wrap_func_expr(mod, expr, wrap_arrays = true)
    @assert expr.head == :function || (expr.head == :(=) &&
                                       expr.args[1] isa Expr &&
                                       expr.args[1].head == :call)

    def = splitdef(expr)

    sig = expr.args[1]
    body = def[:body]

    fname = def[:name]
    args = get(def, :args, [])
    kwargs = get(def, :kwargs, [])

    if fname isa Expr && fname.head == :(::) && length(fname.args) > 1
        self = fname.args[1]
    else
        self = :nothing # LOL -- in this case the argument named nothing is passed nothing
    end
    impl_name = Symbol(fname,"_", hash(string(args)*string(kwargs)))

    function kwargname(kwarg)
        if kwarg isa Expr && kwarg.head == :kw
            kwarg.args[1]
        elseif kwarg isa Expr && kwarg.head == :(...)
            kwarg.args[1]
        else
            kwarg
        end
    end

    function argname(arg)
        if arg isa Expr && (arg.head == :(::) || arg.head == :(...))
            arg.args[1]
        elseif arg isa Expr
            error("$arg not supported as an argument")
        else
            arg
        end
    end

    names = vcat(argname.(args), kwargname.(kwargs))

    function type_options(arg)
        # for every argument find the types that
        # should be allowed as argument. These are:
        #
        # (1) T    (2) wrapper_type(T)    (3) Symbolic{T}
        #
        # However later while emitting methods we omit the one
        # method where all arguments are (1) since those are
        # expected to be defined outside Symbolics
        if arg isa Expr && arg.head == :(::)
            T = Base.eval(mod, arg.args[2])
            Ts = has_symwrapper(T) ? (T, :(Symbolics.SymbolicUtils.Symbolic{<:$T}), wrapper_type(T)) :
                                (T,:(Symbolics.SymbolicUtils.Symbolic{<:$T}))
            if T <: AbstractArray && wrap_arrays
                eT = eltype(T)
                if eT == Any
                    eT = Real
                end
                _arr_type_fn = if hasmethod(ndims, Tuple{Type{T}})
                    (elT) -> :(AbstractArray{T, $(ndims(T))} where {T <: $elT})
                else
                    (elT) -> :(AbstractArray{T} where {T <: $elT})
                end
                if has_symwrapper(eT)
                    Ts = (Ts..., _arr_type_fn(:(Symbolics.SymbolicUtils.Symbolic{<:$eT})), 
                    _arr_type_fn(wrapper_type(eT)))
                else
                    Ts = (Ts..., _arr_type_fn(:(Symbolics.SymbolicUtils.Symbolic{<:$eT})))
                end
            end
            Ts
        elseif arg isa Expr && arg.head == :(...)
            Ts = type_options(arg.args[1])
            map(x->Vararg{x},Ts)
        else
            (Any,)
        end
    end

    types = map(type_options, args)

    impl = :(function $impl_name($self, $(names...))
        $body
    end)
    # TODO: maybe don't drop first lol
    methods = map(Iterators.drop(Iterators.product(types...), 1)) do Ts
        method_args = map(names, Ts) do n, T
            :($n::$T)
        end

        fbody = :(if any($iswrapped, ($(names...),))
                      $wrap($impl_name($self, $([:($unwrap($arg)) for arg in names]...)))
                  else
                      $impl_name($self, $(names...))
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
    end |> esc
end

macro wrapped(expr, wrap_arrays = true)
    wrap_func_expr(__module__, expr, wrap_arrays)
end
