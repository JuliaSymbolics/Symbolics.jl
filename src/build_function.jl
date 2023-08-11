using Base.Threads
using SymbolicUtils.Code: LazyState

abstract type BuildTargets end
struct JuliaTarget <: BuildTargets end
struct StanTarget <: BuildTargets end
struct CTarget <: BuildTargets end
struct MATLABTarget <: BuildTargets end

abstract type ParallelForm end
struct SerialForm <: ParallelForm end

"""
    ShardedForm{multithread}(cutoff, ncalls)

Split a long array construction into nested functions where each function calls
`ncalls` other functions, and the leaf functions populate at most `cutoff` number
of items in the array. If `multithread` is true, uses threading.
"""
struct ShardedForm{multithreaded} <: ParallelForm
    cutoff::Union{Nothing,Int}
    ncalls::Int
end

ShardedForm(cutoff, ncalls) = ShardedForm{false}(cutoff, ncalls)
ShardedForm() = ShardedForm(80, 4)

const MultithreadedForm = ShardedForm{true}

MultithreadedForm() = MultithreadedForm(nothing, 2*nthreads())

function throw_missing_specialization(n)
    throw(ArgumentError("Missing specialization for $n arguments. Check `iip_config`."))
end

"""
`build_function`

Generates a numerically-usable function from a Symbolics `Num`.

```julia
build_function(ex, args...;
               expression = Val{true},
               target = JuliaTarget(),
               parallel=nothing,
               kwargs...)
```

Arguments:

- `ex`: The `Num` to compile
- `args`: The arguments of the function
- `expression`: Whether to generate code or whether to generate the compiled form.
  By default, `expression = Val{true}`, which means that the code for the
  function is returned. If `Val{false}`, then the returned value is compiled.

Keyword Arguments:

- `target`: The output target of the compilation process. Possible options are:
    - `JuliaTarget`: Generates a Julia function
    - `CTarget`: Generates a C function
    - `StanTarget`: Generates a function for compiling with the Stan probabilistic
      programming language
    - `MATLABTarget`: Generates an anonymous function for use in MATLAB and Octave
      environments
- `parallel`: The kind of parallelism to use in the generated function. Defaults
  to `SerialForm()`, i.e. no parallelism, if `ex` is a single expression or an
  array containing <= 1500 non-zero expressions. If `ex` is an array of > 1500
  non-zero expressions, then `ShardedForm(80, 4)` is used. See below for more on
  `ShardedForm`.
  Note that the parallel forms are not exported and thus need to be chosen like
  `Symbolics.SerialForm()`.
  The choices are:
  - `SerialForm()`: Serial execution.
  - `ShardedForm(cutoff, ncalls)`: splits the output function into sub-functions
     which contain at most `cutoff` number of output `rhss`. These sub-functions
     are called by the top-level function that _build_function returns.
     This helps in reducing the compile time of the generated function.
  - `MultithreadedForm()`: Multithreaded execution with a static split, evenly
    splitting the number of expressions per thread.
- `fname`: Used by some targets for the name of the function in the target space.

Note that not all build targets support the full compilation interface. Check the
individual target documentation for details.
"""
function build_function(args...;target = JuliaTarget(),kwargs...)
  _build_function(target,args...;kwargs...)
end

# Scalar output

unwrap_nometa(x) = unwrap(x)
unwrap_nometa(x::CallWithMetadata) = unwrap(x.f)
function destructure_arg(arg::Union{AbstractArray, Tuple,NamedTuple}, inbounds, name)
    if !(arg isa Arr)
        DestructuredArgs(map(unwrap_nometa, arg), name, inbounds=inbounds, create_bindings=false)
    else
        unwrap_nometa(arg)
    end
end
destructure_arg(arg, _, _) = unwrap_nometa(arg)

function _build_function(target::JuliaTarget, op, args...;
                         conv = toexpr,
                         expression = Val{true},
                         expression_module = @__MODULE__(),
                         checkbounds = false,
                         states = LazyState(),
                         linenumbers = true,
                         wrap_code = nothing,
                         cse = false, kwargs...)
  dargs = map((x) -> destructure_arg(x[2], !checkbounds, Symbol("ˍ₋arg$(x[1])")), enumerate([args...]))
    expr = if cse
        fun = Func(dargs, [], Code.cse(unwrap(op)))
        (wrap_code !== nothing) && (fun = wrap_code(fun))
        toexpr(fun, states)
    else
        fun = Func(dargs, [], op)
        (wrap_code !== nothing) && (fun = wrap_code(fun))
        toexpr(fun, states)
    end

    if expression == Val{true}
        expr
    else
        _build_and_inject_function(expression_module, expr)
    end
end

SymbolicUtils.Code.get_symbolify(x::Arr) = SymbolicUtils.Code.get_symbolify(unwrap(x))

function _build_function(target::JuliaTarget, op::Union{Arr, ArrayOp}, args...;
                         conv = toexpr,
                         expression = Val{true},
                         expression_module = @__MODULE__(),
                         checkbounds = false,
                         states = LazyState(),
                         linenumbers = true,
                         cse = false, kwargs...)

    dargs = map((x) -> destructure_arg(x[2], !checkbounds,
                                  Symbol("ˍ₋arg$(x[1])")), enumerate([args...]))

    expr = if cse
        toexpr(Func(dargs, [], Code.cse(unwrap(op))), states)
    else
        toexpr(Func(dargs, [], op), states)
    end

    outsym = Symbol("ˍ₋out")
    body = inplace_expr(unwrap(op), outsym)
    oop_expr = toexpr(Func([outsym, dargs...], [], body), states)

    N = length(shape(op))
    op = unwrap(op)
    if op isa ArrayOp && istree(op.term)
        op_body = op.term
    else
        op_body = :(let $outsym = zeros(Float64, map(length, ($(shape(op)...),)))
                   $body
              $outsym
          end) |> LiteralExpr
    end
    ip_expr = toexpr(Func(dargs, [], op_body), states)
    if expression == Val{true}
        oop_expr, ip_expr
    else
        _build_and_inject_function(expression_module, oop_expr),
        _build_and_inject_function(expression_module, ip_expr)
    end
end

function _build_and_inject_function(mod::Module, ex)
    if ex.head == :function && ex.args[1].head == :tuple
        ex.args[1] = Expr(:call, :($mod.$(gensym())), ex.args[1].args...)
    elseif ex.head == :(->)
        return _build_and_inject_function(mod, Expr(:function, ex.args...))
    end
    # XXX: Workaround to specify the module as both the cache module AND context module.
    # Currently, the @RuntimeGeneratedFunction macro only sets the context module.
    module_tag = getproperty(mod, RuntimeGeneratedFunctions._tagname)
    RuntimeGeneratedFunctions.RuntimeGeneratedFunction(module_tag, module_tag, ex; opaque_closures=false)
end

toexpr(n::Num, st) = toexpr(value(n), st)

function fill_array_with_zero!(x::AbstractArray)
    if eltype(x) <: AbstractArray
        foreach(fill_array_with_zero!, x)
    else
        fill!(x, false)
    end
    return x
end

"""
Build function target: `JuliaTarget`

```julia
function _build_function(target::JuliaTarget, rhss, args...;
                         conv = toexpr, expression = Val{true},
                         checkbounds = false,
                         linenumbers = false,
                         headerfun = addheader, outputidxs=nothing,
                         convert_oop = true, force_SA = false,
                         skipzeros = outputidxs===nothing,
                         fillzeros = skipzeros && !(typeof(rhss)<:SparseMatrixCSC),
                         parallel=SerialForm(), kwargs...)
```

Generates a Julia function which can then be utilized for further evaluations.
If expression=Val{false}, the return is a Julia function which utilizes
RuntimeGeneratedFunctions.jl to be free of world-age issues.

If the `rhss` is a scalar, the generated function is a function
with a scalar output. Otherwise, if it's an `AbstractArray`, the output
is two functions, one for out-of-place AbstractArray output and a second which
is a mutating function. The outputted functions match the given argument order,
i.e., f(u,p,args...) for the out-of-place and scalar functions and
`f!(du,u,p,args..)` for the in-place version.

Special Keyword Arguments:

- `parallel`: The kind of parallelism to use in the generated function. Defaults
  to `SerialForm()`, i.e. no parallelism. Note that the parallel forms are not
  exported and thus need to be chosen like `Symbolics.SerialForm()`.
  The choices are:
  - `SerialForm()`: Serial execution.
  - `ShardedForm(cutoff, ncalls)`: splits the output function into sub-functions
     which contain at most `cutoff` number of output `rhss`. These sub-functions
     are called by the top-level function that _build_function returns.
  - `MultithreadedForm()`: Multithreaded execution with a static split, evenly
    splitting the number of expressions per thread.
- `conv`: The conversion function of symbolic types to Expr. By default, this uses
  the `toexpr` function.
- `checkbounds`: For whether to enable bounds checking inside the generated
  function. Defaults to false, meaning that `@inbounds` is applied.
- `linenumbers`: Determines whether the generated function expression retains
  the line numbers. Defaults to true.
- `convert_oop`: Determines whether the OOP version should try to convert
  the output to match the type of the first input. This is useful for
  cases like LabelledArrays or other array types that carry extra
  information. Defaults to true.
- `force_SA`: Forces the output of the OOP version to be a StaticArray.
  Defaults to `false`, and outputs a static array when the first argument
  is a static array.
- `skipzeros`: Whether to skip filling zeros in the in-place version if the
  filling function is 0.
- `fillzeros`: Whether to perform `fill(out,0)` before the calculations to ensure
  safety with `skipzeros`.
"""
function _build_function(target::JuliaTarget, rhss::AbstractArray, args...;
                       expression = Val{true},
                       expression_module = @__MODULE__(),
                       checkbounds = false,
                       postprocess_fbody=ex -> ex,
                       linenumbers = false,
                       outputidxs=nothing,
                       skipzeros = false,
                       force_SA = false,
                       wrap_code = (nothing, nothing),
                       fillzeros = skipzeros && !(rhss isa SparseMatrixCSC),
                       states = LazyState(),
                       iip_config = (true, true),
                       parallel=nothing, cse = false, kwargs...)

    # We cannot switch to ShardedForm because it deadlocks with
    # RuntimeGeneratedFunctions
    dargs = map((x) -> destructure_arg(x[2], !checkbounds,
                                  Symbol("ˍ₋arg$(x[1])")), enumerate([args...]))
    i = findfirst(x->x isa DestructuredArgs, dargs)
    similarto = force_SA ? SArray : i === nothing ? Array : dargs[i].name

    oop, iip = iip_config
    oop_body = if oop
        postprocess_fbody(make_array(parallel, dargs, rhss, similarto, cse))
    else
        term(throw_missing_specialization, length(dargs))
    end
    oop_expr = Func(dargs, [], oop_body)

    if !isnothing(wrap_code[1])
        oop_expr = wrap_code[1](oop_expr)
    end

    out = Sym{Any}(:ˍ₋out)
    ip_body = if iip
        postprocess_fbody(set_array(parallel,
                                    dargs,
                                    out,
                                    outputidxs,
                                    rhss,
                                    checkbounds,
                                    skipzeros,
                                    cse,))
    else
        term(throw_missing_specialization, length(dargs) + 1)
    end
    ip_expr = Func([out, dargs...], [], ip_body)

    if !isnothing(wrap_code[2])
        ip_expr = wrap_code[2](ip_expr)
    end

    if expression == Val{true}
        return toexpr(oop_expr, states), toexpr(ip_expr, states)
    else
        return _build_and_inject_function(expression_module, toexpr(oop_expr, states)),
        _build_and_inject_function(expression_module, toexpr(ip_expr, states))
    end
end

_nnz(x::AbstractArray) = length(x)
_nnz(x::AbstractSparseArray) = nnz(x)
_nnz(x::Union{Base.ReshapedArray, LinearAlgebra.Transpose}) = _nnz(parent(x))

function make_array(s, dargs, arr, similarto, cse)
    s !== nothing && Base.@warn("Parallel form of $(typeof(s)) not implemented")
    _make_array(arr, similarto, cse)
end

function make_array(s::SerialForm, dargs, arr, similarto, cse)
    _make_array(arr, similarto, cse)
end

function make_array(s::ShardedForm, closed_args, arr, similarto, cse)
    if arr isa AbstractSparseArray

        return LiteralExpr(quote
                               $SparseMatrixCSC($(arr.m),
                                               $(arr.n),
                                               copy($(arr.colptr)),
                                               copy($(arr.rowval)),
                                               $(make_array(s,
                                                            closed_args,
                                                            arr.nzval,
                                                            Vector,cse)))
                           end)
    end
    per_task = ceil(Int, length(arr) / s.ncalls)
    slices = collect(Iterators.partition(arr, per_task))
    arrays = map(slices) do slice
        Func(closed_args, [], _make_array(slice, similarto, cse)), closed_args
    end
    SpawnFetch{typeof(s)}(first.(arrays), last.(arrays), vcat)
end

struct Funcall{F, T}
    f::F
    args::T
end

(f::Funcall)() = f.f(f.args...)

function toexpr(p::SpawnFetch{MultithreadedForm}, st)
    args = isnothing(p.args) ?
              Iterators.repeated((), length(p.exprs)) : p.args
    spawns = map(p.exprs, args) do thunk, a
        ex = :($Funcall($(drop_expr(@RuntimeGeneratedFunction(@__MODULE__, toexpr(thunk, st), false))),
                       ($(toexpr.(a, (st,))...),)))
        quote
            let
                task = Base.Task($ex)
                task.sticky = false
                Base.schedule(task)
                task
            end
        end
    end
    quote
        $(toexpr(p.combine, st))(map(fetch, ($(spawns...),))...)
    end
end

function toexpr(p::SpawnFetch{ShardedForm{false}}, st)
    args = isnothing(p.args) ?
              Iterators.repeated((), length(p.exprs)) : p.args
    spawns = map(p.exprs, args) do thunk, a
        :($(drop_expr(@RuntimeGeneratedFunction(@__MODULE__, toexpr(thunk, st), false)))($(toexpr.(a, (st,))...),))
    end
    quote
        $(toexpr(p.combine, st))($(spawns...))
    end
end

function nzmap(f, x::Union{Base.ReshapedArray, LinearAlgebra.Transpose})
    Setfield.@set x.parent = nzmap(f, x.parent)
end

function nzmap(f, x::SubArray)
    unview = copy(x)
    if unview isa Union{SparseMatrixCSC, SparseVector}
        n = nnz(unview)
        if n != length(unview.nzval)
            resize!(unview.nzval, n)
            resize!(unview.rowval, n)
        end
    end
    nzmap(f, unview)
end

function nzmap(f, x::AbstractSparseArray)
    Setfield.@set x.nzval = nzmap(f, x.nzval)
end
nzmap(f, x) = map(f, x)

_issparse(x::AbstractArray) = issparse(x)
_issparse(x::Union{SubArray, Base.ReshapedArray, LinearAlgebra.Transpose}) = _issparse(parent(x))

function _make_sparse_array(arr, similarto, cse)
    if arr isa Union{SubArray, Base.ReshapedArray, LinearAlgebra.Transpose}
        LiteralExpr(quote
            $Setfield.@set $(nzmap(x->true, arr)).parent =
                $(_make_array(parent(arr), typeof(parent(arr)), cse))
            end)
    else
        LiteralExpr(quote
                        let __reference = copy($(nzmap(x->true, arr)))
                            $Setfield.@set __reference.nzval =
                            $(_make_array(arr.nzval, Vector{symtype(eltype(arr))}, cse))
                        end
                    end)
    end
end

function _make_array(rhss::AbstractArray, similarto, cse)
    arr = nzmap(x->_make_array(x, similarto, cse), rhss)
    if _issparse(arr)
        _make_sparse_array(arr, similarto, cse)
    elseif cse
        Code.cse(MakeArray(unwrap.(arr), similarto))
    else
        MakeArray(arr, similarto)
    end
end

_make_array(x, similarto, cse) = x

## In-place version

function set_array(p, closed_vars, args...)
    p !== nothing && Base.@warn("Parallel form of $(typeof(p)) not implemented")
    _set_array(args...)
end

function set_array(s::SerialForm, closed_vars, args...)
    _set_array(args...)
end

function recursive_split(leaf_f, s, out, args, outputidxs, xs)
    cutoff = isnothing(s.cutoff) ? ceil(Int, length(xs) / (2*s.ncalls)) : s.cutoff
    if length(xs) <= cutoff
        return leaf_f(outputidxs, xs)
    else
        per_part = ceil(Int, length(xs) / s.ncalls)
        slices = collect(Iterators.partition(zip(outputidxs, xs), per_part))
        fs = map(slices) do slice
            recursive_split(leaf_f, s, out, args, first.(slice), last.(slice))
        end
        return Func(args, [],
                    SpawnFetch{typeof(s)}(fs, [args for f in fs],
                                          (@inline noop(x...) = nothing)),
                    [])
    end
end

function set_array(s::ShardedForm, closed_args, out, outputidxs, rhss, checkbounds, skipzeros, cse)
    if rhss isa AbstractSparseArray
        return set_array(s,
                         closed_args,
                         LiteralExpr(:($out.nzval)),
                         nothing,
                         rhss.nzval,
                         checkbounds,
                         skipzeros,
                         cse)
    end

    outvar = !(out isa Sym) ? gensym("out") : out

    if outputidxs === nothing
        outputidxs = collect(eachindex(rhss))
    end
    all_args = [outvar, closed_args...]
    ex = recursive_split(s, outvar, all_args, outputidxs, rhss) do idxs, xs
        Func(all_args, [],
             _set_array(outvar, idxs, xs, checkbounds, skipzeros, cse),
             [])
    end.body

    return out isa Sym ? ex : LiteralExpr(quote
        $outvar = $out
        $ex
    end)
end

function _set_array(out, outputidxs, rhss::AbstractSparseArray, checkbounds, skipzeros, cse)
    _set_array(LiteralExpr(:($out.nzval)), nothing, rhss.nzval, checkbounds, skipzeros, cse)
end

function _set_array(out, outputidxs, rhss::AbstractArray, checkbounds, skipzeros, cse)
    if outputidxs === nothing
        outputidxs = collect(eachindex(rhss))
    end
    # sometimes outputidxs is a Tuple
    ii = findall(i->!(rhss[i] isa AbstractArray) && !(skipzeros && _iszero(rhss[i])), eachindex(outputidxs))
    jj = findall(i->rhss[i] isa AbstractArray, eachindex(outputidxs))
    exprs = []
    setterexpr = SetArray(!checkbounds,
                          out,
                          [AtIndex(outputidxs[i],
                                   rhss[i])
                           for i in ii])
    cse ? push!(exprs, Code.cse(setterexpr)) : push!(exprs, setterexpr)
    for j in jj
        push!(exprs, _set_array(LiteralExpr(:($out[$j])), nothing, rhss[j], checkbounds, skipzeros, cse))
    end
    LiteralExpr(quote
                    $(exprs...)
                end)
end

_set_array(out, outputidxs, rhs, checkbounds, skipzeros, cse) = rhs


function vars_to_pairs(name,vs::Union{Tuple, AbstractArray}, symsdict=Dict())
    vs_names = tosymbol.(vs)
    for (v,k) in zip(vs_names, vs)
        symsdict[k] = Sym{symtype(k)}(v)
    end
    exs = [:($name[$i]) for (i, u) ∈ enumerate(vs)]
    vs_names,exs
end
function vars_to_pairs(name,vs, symsdict)
    symsdict[vs] = Sym{symtype(vs)}(tosymbol(vs))
    [tosymbol(vs)], [name]
end

get_varnumber(varop, vars::Vector) =  findfirst(x->isequal(x,varop),vars)
get_varnumber(varop, var) =  isequal(var,varop) ? 0 : nothing

function buildvarnumbercache(args...)
    varnumsdict = Pair[]
    for (argi,arg) in enumerate(args)
        if isa(arg,AbstractArray)
            for (eli,el) in enumerate(arg)
                push!(varnumsdict, el=>(argi,eli))
            end
        else
            push!(varnumsdict ,arg=>(argi,0))
        end
    end
    return Dict(varnumsdict)
end

function numbered_expr(O::Symbolic,varnumbercache,args...;varordering = args[1],offset = 0,
                       states = LazyState(),
                       lhsname=:du,rhsnames=[Symbol("MTK$i") for i in 1:length(args)])
    O = value(O)
    if (issym(O) || issym(operation(O))) || (istree(O) && operation(O) == getindex)
        (j,i) = get(varnumbercache, O, (nothing, nothing))
        if !isnothing(j)
            return i==0 ? :($(rhsnames[j])) : :($(rhsnames[j])[$(i+offset)])
        end
    end
    if istree(O)
        if operation(O) === getindex
            args = arguments(O)
            Expr(:ref, toexpr(args[1], states), toexpr.(args[2:end] .+ offset, (states,))...)
        else
            Expr(:call, Symbol(operation(O)), (numbered_expr(x,varnumbercache,args...;offset=offset,lhsname=lhsname,
                                                             rhsnames=rhsnames,varordering=varordering) for x in arguments(O))...)
        end
    elseif issym(O)
        tosymbol(O, escape=false)
    else
        O
    end
end

function numbered_expr(de::Equation,varnumbercache,args...;varordering = args[1],
                       lhsname=:du,rhsnames=[Symbol("MTK$i") for i in 1:length(args)],offset=0)

    varordering = value.(args[1])
    var = var_from_nested_derivative(de.lhs)[1]
    i = findfirst(x->isequal(tosymbol(issym(x) ? x : operation(x), escape=false), tosymbol(var, escape=false)),varordering)
    :($lhsname[$(i+offset)] = $(numbered_expr(de.rhs,varnumbercache,args...;offset=offset,
                                              varordering = varordering,
                                              lhsname = lhsname,
                                              rhsnames = rhsnames)))
end
numbered_expr(c,args...;kwargs...) = c
numbered_expr(c::Num,args...;kwargs...) = error("Num found")


# Replace certain multiplication and power expressions so they form valid C code
# Extra factors of 1 are hopefully eliminated by the C compiler
function coperators(expr)
    expr isa Expr || return expr
    for e in expr.args
        if e isa Expr
            coperators(e)
        end
    end
    for i in eachindex(expr.args)
        if expr.args[i] isa Rational
            expr.args[i] = float(expr.args[i]) # Evaluate rational numbers to floating-point
        end
    end
    # Introduce another factor 1 to prevent contraction of terms like "5 * t" to "5t" (not valid C code)
    if expr.head==:call && expr.args[1]==:* && length(expr.args)==3 && isa(expr.args[2], Real) && isa(expr.args[3], Symbol)
        push!(expr.args, 1)
    # Power operator does not exist in C, replace by multiplication or "pow"
    elseif expr.head==:call && expr.args[1]==:^
        @assert length(expr.args)==3 "Don't know how to handle ^ operation with <> 2 arguments"
        x = expr.args[2]
        n = expr.args[3]
        empty!(expr.args)
        # Replace by multiplication/division if
        #   x is a symbol and  n is a small integer
        #   x is a more complex expression and n is ±1
        #   n is exactly 0
        if (isa(n,Integer) && ((isa(x, Symbol) && abs(n) <= 3) || abs(n) <= 1)) || n==0
            if n >= 0
                append!(expr.args, [:*, fill(x, n)...])
                # fill up with factor 1 so this expr can still be a multiplication
                while length(expr.args) < 3
                    push!(expr.args, 1)
                end
            else # inverse of the above
                if n==-1
                    term = x
                else
                    term = :( ($(x)) ^ ($(-n)))
                    coperators(term)
                end
                append!(expr.args, [:/, 1., term])
            end
        #... otherwise use "pow" function
        else
            append!(expr.args, [:pow, x, n])
        end
    end
    expr
end


"""
Build function target: `CTarget`

```julia
function _build_function(target::CTarget, eqs::Array{<:Equation}, args...;
                         conv = toexpr, expression = Val{true},
                         fname = :diffeqf,
                         lhsname=:du,rhsnames=[Symbol("RHS\$i") for i in 1:length(args)],
                         libpath=tempname(),compiler=:gcc)
```

This builds an in-place C function. Only works on arrays of equations. If
`expression == Val{false}`, then this builds a function in C, compiles it,
and returns a lambda to that compiled function. These special keyword arguments
control the compilation:

- libpath: the path to store the binary. Defaults to a temporary path.
- compiler: which C compiler to use. Defaults to :gcc, which is currently the
  only available option.
"""
function _build_function(target::CTarget, eqs::Array{<:Equation}, args...;
                         conv = toexpr, expression = Val{true},
                         fname = :diffeqf,
                         lhsname=:du,rhsnames=[Symbol("RHS$i") for i in 1:length(args)],
                         libpath=tempname(),compiler=:gcc)

    @warn "build_function(::Array{<:Equation}...) is deprecated. Use build_function(::AbstractArray...) instead."

    varnumbercache = buildvarnumbercache(args...)
    differential_equation = string(join([numbered_expr(eq,varnumbercache,args...,lhsname=lhsname,
                                  rhsnames=rhsnames,offset=-1) for
                                  (i, eq) ∈ enumerate(eqs)],";\n  "),";")

    argstrs = join(vcat("double* $(lhsname)",[typeof(args[i])<:AbstractArray ? "double* $(rhsnames[i])" : "double $(rhsnames[i])" for i in 1:length(args)]),", ")
    ex = """
    void $fname($(argstrs...)) {
      $differential_equation
    }
    """

    if expression == Val{true}
        return ex
    else
        @assert compiler == :gcc
        ex = build_function(eqs,args...;target=Symbolics.CTarget())
        open(`gcc -fPIC -O3 -msse3 -xc -shared -o $(libpath * "." * Libdl.dlext) -`, "w") do f
            print(f, ex)
        end
        drop_expr(@RuntimeGeneratedFunction(@__MODULE__,
                                            :((du::Array{Float64},u::Array{Float64},p::Array{Float64},t::Float64)
                                              -> ccall(("diffeqf", $libpath),
                                                       Cvoid, (Ptr{Float64},
                                                               Ptr{Float64},
                                                               Ptr{Float64},
                                                               Float64), du, u,
                                                               p, t)), false))
    end
end


"""
Build function target: `CTarget`

```julia
function _build_function(target::CTarget, ex::AbstractArray, args...;
                         columnmajor = true,
                         conv        = toexpr,
                         expression  = Val{true},
                         fname       = :diffeqf,
                         lhsname     = :du,
                         rhsnames    = [Symbol("RHS\$i") for i in 1:length(args)],
                         libpath     = tempname(),
                         compiler    = :gcc)
```

This builds an in-place C function. Only works on expressions. If
`expression == Val{false}`, then this builds a function in C, compiles it,
and returns a lambda to that compiled function. These special keyword arguments
control the compilation:

- libpath: the path to store the binary. Defaults to a temporary path.
- compiler: which C compiler to use. Defaults to :gcc, which is currently the
  only available option.
"""
function _build_function(target::CTarget, ex::AbstractArray, args...;
                         columnmajor = true,
                         conv        = toexpr,
                         expression  = Val{true},
                         fname       = :diffeqf,
                         lhsname     = :du,
                         rhsnames    = [Symbol("RHS$i") for i in 1:length(args)],
                         libpath     = tempname(),
                         compiler    = :gcc)

    if !columnmajor
        return _build_function(target, hcat([row for row ∈ eachrow(ex)]...), args...;
                               columnmajor = true,
                               conv        = conv,
                               fname       = fname,
                               lhsname     = lhsname,
                               rhsnames    = rhsnames,
                               libpath     = libpath,
                               compiler    = compiler)
    end


    varnumbercache = buildvarnumbercache(args...)
    equations = Vector{String}()
    for col ∈ 1:size(ex,2)
        for row ∈ 1:size(ex,1)
            lhs = string(lhsname, "[", (col-1) * size(ex,1) + row-1, "]")
            rhs = numbered_expr(value(ex[row, col]), varnumbercache, args...;
                                lhsname  = lhsname,
                                rhsnames = rhsnames,
                                offset   = -1) |> coperators |> string  # Filter through coperators to produce valid C code in more cases
            push!(equations, string(lhs, " = ", rhs, ";"))
        end
    end

    argstrs = join(vcat("double* $(lhsname)",[typeof(args[i])<:AbstractArray ? "const double* $(rhsnames[i])" : "const double $(rhsnames[i])" for i in 1:length(args)]),", ")

    ccode = """
    #include <math.h>
    void $fname($(argstrs...)) {$([string("\n  ", eqn) for eqn ∈ equations]...)\n}
    """

    if expression == Val{true}
        return ccode
    else
        @assert compiler == :gcc
        open(`gcc -fPIC -O3 -msse3 -xc -shared -o $(libpath * "." * Libdl.dlext) -`, "w") do f
            print(f, ccode)
        end
        drop_expr(@RuntimeGeneratedFunction(@__MODULE__,
                                            :((du::Array{Float64},u::Array{Float64},p::Array{Float64},t::Float64)
                                              -> ccall(("diffeqf", $libpath),
                                                       Cvoid, (Ptr{Float64},
                                                               Ptr{Float64},
                                                               Ptr{Float64},
                                                               Float64), du, u,
                                                       p, t)), false))
    end

end
_build_function(target::CTarget, ex::Num, args...; kwargs...) = _build_function(target, [ex], args...; kwargs...)


"""
Build function target: `StanTarget`

```julia
function _build_function(target::StanTarget, eqs::Array{<:Equation}, vs, ps, iv;
                         conv = toexpr, expression = Val{true},
                         fname = :diffeqf, lhsname=:internal_var___du,
                         rhsnames=[:internal_var___u,:internal_var___p,:internal_var___t])
```

This builds an in-place Stan function compatible with the Stan differential equation solvers.
Unlike other build targets, this one requires (vs, ps, iv) as the function arguments.
Only allowed on arrays of equations.
"""
function _build_function(target::StanTarget, eqs::Array{<:Equation}, vs, ps, iv;
                         conv = toexpr, expression = Val{true},
                         fname = :diffeqf, lhsname=:internal_var___du,
                         rhsnames=[:internal_var___u,:internal_var___p,:internal_var___t])

    @warn "build_function(::Array{<:Equation}...) is deprecated. Use build_function(::AbstractArray...) instead."
    @assert expression == Val{true}

    varnumbercache = buildvarnumbercache(vs,ps...)
    par_str = join(["real $(rhsnames[2])_$i" for i in 1:length(ps)], ", ")
    rhsnames_mod = [:internal_var___u, [Symbol("$(rhsnames[2])_$i") for i in 1:length(ps)]..., :internal_var___t]
    differential_equation = string(join([numbered_expr(eq,varnumbercache,vs,ps,lhsname=lhsname,
                                   rhsnames=rhsnames_mod) for
                                   (i, eq) ∈ enumerate(eqs)],";\n  "),";")


    """
    vector $fname(real $(conv(iv)),vector $(rhsnames[1]), $par_str) {
      vector[$(length(eqs))] $lhsname;
      $differential_equation
      return $lhsname;
    }
    """
end

"""
Build function target: `StanTarget`

```julia
function _build_function(target::StanTarget, ex::AbstractArray, vs, ps, iv;
                         columnmajor = true,
                         conv        = toexpr,
                         expression  = Val{true},
                         fname       = :diffeqf, lhsname=:internal_var___du,
                         rhsnames    =  [:internal_var___u,:internal_var___p,:internal_var___t])
```

This builds an in-place Stan function compatible with the Stan differential equation solvers.
Unlike other build targets, this one requires (vs, ps, iv) as the function arguments.
Only allowed on expressions, and arrays of expressions.
"""
function _build_function(target::StanTarget, ex::AbstractArray, vs, ps, iv;
                         columnmajor = true,
                         conv        = toexpr,
                         expression  = Val{true},
                         fname       = :diffeqf, lhsname=:internal_var___du,
                         rhsnames    =  [:internal_var___u,:internal_var___p,:internal_var___t])

    @assert expression == Val{true}

    if !columnmajor
        return _build_function(target, hcat([row for row ∈ eachrow(ex)]...), vs, ps, iv;
                            columnmajor = true,
                            conv        = conv,
                            expression  = expression,
                            fname       = fname,
                            lhsname     = lhsname,
                            rhsnames    = rhsnames)
    end

    varnumbercache = buildvarnumbercache(vs,ps,iv)
    equations = Vector{String}()
    for col ∈ 1:size(ex,2)
        for row ∈ 1:size(ex,1)
            lhs = string(lhsname, "[", (col-1) * size(ex,1) + row, "]")
            rhs = numbered_expr(value(ex[row, col]), varnumbercache, vs, ps, iv;
                                lhsname  = lhsname,
                                rhsnames = rhsnames,
                                offset   = 0) |> string
            push!(equations, string(lhs, " = ", rhs, ";"))
        end
    end

    """
    vector $fname(real $(conv(iv)),vector $(rhsnames[1]),vector $(rhsnames[2])) {
      vector[$(length(equations))] $lhsname;
    $([eqn == equations[end] ? string("  ", eqn) : string("  ", eqn, "\n") for eqn ∈ equations]...)
      return $lhsname;
    }
    """
end
_build_function(target::StanTarget, ex::Num, vs, ps, iv; kwargs...) = _build_function(target, [ex], vs, ps, iv; kwargs...)

"""
Build function target: `MATLABTarget`

```julia
function _build_function(target::MATLABTarget, eqs::Array{<:Equation}, args...;
                         conv = toexpr, expression = Val{true},
                         lhsname=:internal_var___du,
                         rhsnames=[:internal_var___u,:internal_var___p,:internal_var___t])
```

This builds an out of place anonymous function @(t,rhsnames[1]) to be used in MATLAB.
Compatible with the MATLAB differential equation solvers. Only allowed on expressions,
and arrays of expressions.
"""
function _build_function(target::MATLABTarget, eqs::Array{<:Equation}, args...;
                         conv = toexpr, expression = Val{true},
                         fname = :diffeqf, lhsname=:internal_var___du,
                         rhsnames=[:internal_var___u,:internal_var___p,:internal_var___t])

    @warn "build_function(::Array{<:Equation}...) is deprecated. Use build_function(::AbstractArray...) instead."
    @assert expression == Val{true}

    varnumbercache = buildvarnumbercache(args...)
    matstr = join([numbered_expr(eq.rhs,varnumbercache,args...,lhsname=lhsname,
                                  rhsnames=rhsnames) for
                                  (i, eq) ∈ enumerate(eqs)],"; ")

    matstr = replace(matstr,"["=>"(")
    matstr = replace(matstr,"]"=>")")
    matstr = "$fname = @($(rhsnames[3]),$(rhsnames[1])) ["*matstr*"];"
    matstr
end

"""
Build function target: `MATLABTarget`

```julia
function _build_function(target::MATLABTarget, ex::AbstractArray, args...;
                         columnmajor = true,
                         conv        = toexpr,
                         expression  = Val{true},
                         fname       = :diffeqf,
                         lhsname     = :internal_var___du,
                         rhsnames    = [:internal_var___u,:internal_var___p,:internal_var___t])
```

This builds an out of place anonymous function @(t,rhsnames[1]) to be used in MATLAB.
Compatible with the MATLAB differential equation solvers. Only allowed on expressions,
and arrays of expressions.
"""
function _build_function(target::MATLABTarget, ex::AbstractArray, args...;
                         columnmajor = true,
                         conv        = toexpr,
                         expression  = Val{true},
                         fname       = :diffeqf,
                         lhsname     = :internal_var___du,
                         rhsnames    = [:internal_var___u,:internal_var___p,:internal_var___t])

    @assert expression == Val{true}

    if !columnmajor
        return _build_function(target, hcat([row for row ∈ eachrow(ex)]...), args...;
                               columnmajor = true,
                               conv        = conv,
                               expression  = expression,
                               fname       = fname,
                               lhsname     = lhsname,
                               rhsnames    = rhsnames)
    end

    varnumbercache = buildvarnumbercache(args...)
    matstr = ""
    for row ∈ 1:size(ex,1)
        row_strings = Vector{String}()
        for col ∈ 1:size(ex,2)
            lhs = string(lhsname, "[", (col-1) * size(ex,1) + row-1, "]")
            rhs = numbered_expr(value(ex[row, col]), varnumbercache, args...;
                                lhsname  = lhsname,
                                rhsnames = rhsnames,
                                offset   = 0) |> string
            push!(row_strings, rhs)
        end
        matstr = matstr * "  " * join(row_strings, ", ") * ";\n"
    end

    matstr = replace(matstr,"["=>"(")
    matstr = replace(matstr,"]"=>")")
    matstr = "$fname = @($(rhsnames[3]),$(rhsnames[1])) [\n"*matstr*"];\n"

    return matstr

end
_build_function(target::MATLABTarget, ex::Num, args...; kwargs...) = _build_function(target, [ex], args...; kwargs...)
