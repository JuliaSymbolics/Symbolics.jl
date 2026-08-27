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

function _get_type_name(x::Union{Symbol, Expr})
    x isa Symbol && return x
    @assert Meta.isexpr(x, :curly)
    return x.args[1]
end

"""
    @symbolic_wrap struct W <: T
        ...
    end

Define `W` as the symbolic wrapper type for the symbolic type `T`, and register the
correspondence between the two with Symbolics.

Symbolic expressions are stored as untyped expression trees (`BasicSymbolic`), but Julia
code dispatches on types: a function written for `Real` will not accept an expression tree.
Symbolics bridges the two by *wrapping* an expression tree in a struct that subtypes the
type the expression stands for. `Num <: Real` is the built-in example — it holds an
expression tree whose `symtype` is `Real` and can therefore be passed anywhere a `Real` is
expected. `@symbolic_wrap` is how a package declares its own such pairing, so that
Symbolics knows to wrap results of that symbolic type in `W` and to unwrap `W` back to the
expression tree at the boundaries.

The macro takes a struct definition whose supertype is the symbolic type being wrapped. It
emits the struct unchanged, plus the trait methods that tie `W` and `T` together:
`Symbolics.has_symwrapper(::Type{<:T})`, `Symbolics.wrapper_type(::Type{<:T}) = W`,
`Symbolics.is_wrapper_type(::Type{<:W})`, `Symbolics.wraps_type(::Type{W}) = T` and
`Symbolics.iswrapped(::W)`. Once those exist, [`Symbolics.wrap`](@ref) maps a value of
symbolic type `T` to a `W`, [`SymbolicUtils.unwrap`](https://symbolicutils.juliasymbolics.org/api/#SymbolicUtils.unwrap) maps it back, [`@wrapped`](@ref)
generates wrapper-accepting methods, and `@register_symbolic` knows `W` counts as symbolic.

Two things are the caller's responsibility. `W` must be constructible from the value it
wraps, because that is how `Symbolics.wrap` builds it, and a `Symbolics.unwrap` method must
be defined for `W` to get the value back out.

The registration is on the *unparameterized* supertype, so
`Symbolics.wrapper_type(T{Int})` also returns `W`; if a parameterized wrapper is wanted
there, add the method by hand.

Only one wrapper may be registered per symbolic type, and `Real` is already taken by `Num`.
`@symbolic_wrap` is therefore for introducing wrappers over *new* symbolic types, not for
replacing the built-in ones.

# Example

```julia
using Symbolics

abstract type AbstractFoo{T} end
struct Foo{T} <: AbstractFoo{T} end

@symbolic_wrap struct FooWrap{T} <: AbstractFoo{T}
    val::Foo{T}
end

Symbolics.unwrap(r::FooWrap) = r.val

wrap(Foo{Int}())            # FooWrap{Int}
unwrap(wrap(Foo{Int}()))    # Foo{Int}
```

With that in place, [`@wrapped`](@ref) can generate methods that accept `FooWrap` wherever
they declare an `AbstractFoo` argument.

See also: [`@wrapped`](@ref), [`Symbolics.wrap`](@ref),
[`SymbolicUtils.unwrap`](https://symbolicutils.juliasymbolics.org/api/#SymbolicUtils.unwrap).
"""
macro symbolic_wrap(expr)
    @assert expr isa Expr && expr.head == :struct
    @assert expr.args[2].head == :(<:)
    sig = expr.args[2]
    T = _get_type_name(sig.args[1])
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

"""
    $(TYPEDSIGNATURES)

Wrap the symbolic or non-symbolic value `x` in the appropriate wrapper type.
"""
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
        # (1) T    (2) wrapper_type(T)    (3) BasicSymbolic
        #
        # However later while emitting methods we omit the one
        # method where all arguments are (1) since those are
        # expected to be defined outside Symbolics
        if arg isa Expr && arg.head == :(::)
            T = Base.eval(mod, arg.args[2])
            if has_symwrapper(T)
                Ts = (T, SymbolicT, wrapper_type(T))
            else
                Ts = (T, SymbolicT)
            end
            if T <: AbstractArray && wrap_arrays
                eT = eltype(T)
                if eT == Any
                    eT = Real
                end
                _arr_type_fn = if hasmethod(ndims, Tuple{Type{T}})
                    (elT) -> AbstractArray{S, ndims(T)} where {S <: elT}
                else
                    (elT) -> AbstractArray{S} where {S <: elT}
                end
                if has_symwrapper(eT)
                    Ts = (Ts..., _arr_type_fn(SymbolicT), _arr_type_fn(wrapper_type(eT)))
                else
                    Ts = (Ts..., _arr_type_fn(SymbolicT))
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

    function is_wrapped_array_eltype(T)
        T <: AbstractArray || return false
        if T isa UnionAll
            _inner = T.body
            while _inner isa UnionAll
                _inner = _inner.body
            end
            return _inner.parameters[1] == T.var && is_wrapper_type(T.var.ub)::Bool ||
                is_wrapper_type(eltype(T))::Bool
        elseif T isa DataType
            return is_wrapper_type(eltype(T))::Bool
        end
        error("Unreachable")
    end
    # TODO: maybe don't drop first lol
    methods = map(Iterators.drop(Iterators.product(types...), 1)) do Ts
        method_args = map(names, Ts) do n, T
            :($n::$T)
        end

        any_wrapper = false
        impl_args = map(enumerate(names)) do (i, name)
            if is_wrapper_type(Ts[i])
                any_wrapper = true
                :($unwrap($name))
            elseif Ts[i] <: AbstractArray && is_wrapped_array_eltype(Ts[i])
                any_wrapper = true
                :($_recursive_unwrap($name))
            else
                name
            end
        end
        implcall = :($impl_name($self, $(impl_args...)))
        if any_wrapper
            implcall = :($wrap($implcall))
        end

        body = Expr(:block)
        for (i, T) in enumerate(Ts)
            if T === BasicSymbolic{VartypeT}
                push!(body.args, :(@assert $symtype($(names[i])) <: $(types[i][1])))
            elseif T <: (AbstractArray{S} where {S <: SymbolicT}) && eltype(types[i][1]) !== Any
                push!(body.args, :(@assert $symtype($(names[i])[1]) <: $(eltype(types[i][1]))))
            end
        end
        push!(body.args, implcall)

        if isempty(kwargs)
            :(function $fname($(method_args...))
                  $body
              end)
        else
            :(function $fname($(method_args...); $(kwargs...))
                  $body
              end)
        end
    end

    quote
        $impl
        $(methods...)
    end |> esc
end

"""
    @wrapped function f(x::T, ...) ... end
    @wrapped f(x::T, ...) = ...
    @wrapped function f(x::T, ...) ... end false

Given one function definition written against concrete types, generate the additional
methods that accept symbolic expressions and wrapper types in those argument positions.

A function written as `f(x::Real)` will not accept a `Num`, because `Num` holds an
expression tree rather than a number, and will not accept a raw `BasicSymbolic` either.
Writing the symbolic methods by hand means one method per combination of "concrete /
symbolic / wrapped" across all arguments, each of which has to unwrap its wrapper
arguments, call the body, and re-wrap the result. `@wrapped` writes that combinatorial
expansion for you.

The body is emitted once under a private name, and one method of `f` is emitted per element
of the product of the type options for each argument. For an argument annotated `::T` the
options are `T` itself, a symbolic expression whose `symtype` is `T`, and — when `T` has a
wrapper registered via [`@symbolic_wrap`](@ref) — that wrapper type. When `T` is an array
type, array-of-symbolic and array-of-wrapper options are added as well. Unannotated
arguments are left as `Any`.

The all-concrete method is deliberately *not* emitted: `@wrapped function f(x::Real)` does
not define `f(::Real)`. That method is assumed to already exist outside Symbolics (often it
is the very Base function being extended), and emitting it would be piracy. Consequently
`f` must have a non-symbolic definition of its own if it is to be called on plain values.

In each generated method, wrapper-typed arguments are unwrapped before the body runs, and
the result is re-wrapped with [`Symbolics.wrap`](@ref) whenever at least one argument was
wrapped — so passing `Num`s in gets a `Num` back out, while passing raw `BasicSymbolic`s in
returns a raw expression. Symbolic arguments additionally get an `@assert` that their
`symtype` matches the declared annotation, which turns a type mismatch into an error at the
call rather than a wrong answer downstream.

Keyword arguments are forwarded as declared and are not expanded over; only positional
arguments participate in the product.

The optional trailing argument (default `true`) controls whether array-typed arguments are
expanded over their symbolic-array and wrapped-array options. Pass `false` to suppress
that, which is useful when the method is meant to see the container itself rather than a
symbolic array.

# Example

Continuing the wrapper from [`@symbolic_wrap`](@ref):

```julia
@wrapped function foo(f::AbstractFoo, x::Real)
    x > 1 ? Foo{Int}() : 0
end

applicable(foo, Foo{Int}(), 2)             # false: the all-concrete method is not emitted
applicable(foo, Foo{Int}(), wrap(2))       # true
applicable(foo, wrap(Foo{Int}()), 2)       # true
applicable(foo, wrap(Foo{Int}()), wrap(2)) # true
```

See also: [`@symbolic_wrap`](@ref), [`Symbolics.wrap`](@ref),
[`SymbolicUtils.unwrap`](https://symbolicutils.juliasymbolics.org/api/#SymbolicUtils.unwrap).
"""
macro wrapped(expr, wrap_arrays = true)
    wrap_func_expr(__module__, expr, wrap_arrays)
end
