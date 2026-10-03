struct FunctionUnimplementedError <: Exception
    type::String
end

function Base.showerror(io::IO, err::FunctionUnimplementedError)
    print(io, "The $(err.type) version for this function is invalid")
end

function get_unimplemented_fn(nargs, type)
    expr = :(function unimplemented()
        throw($FunctionUnimplementedError($type))
     end)

    for i in 1:nargs
        push!(expr.args[1].args, Symbol(:x, i))
    end
    return expr
end

const IIP_OUTSYM = only(@syms $DEFAULT_OUTSYM::Any)
const IIP_ALLOCATOR = SU.Term{VartypeT}(
    Returns, SArgsT((IIP_OUTSYM,));
    type = FnType{Tuple, Any, Returns{Any}}, shape = SU.ShapeVecT()
)

function canonicalize_args(args::Vector, inbounds::Bool)
    return map(enumerate(args)) do (i, arg)
        if arg isa Arr
            unwrap(arg)
        elseif arg isa AbstractArray
            DestructuredArgs(map(unwrap, arg), default_arg_name(i); inbounds, create_bindings = false)
        elseif arg isa Union{Tuple, NamedTuple}
            DestructuredArgs(map(unwrap, collect(arg)), default_arg_name(i); inbounds, create_bindings = false)
        else
            unwrap(arg)
        end
    end
end

"""
    CodegenFunctionOptions(; kwargs...)

Bundle of options controlling code generation in [`codegen_function`](@ref).

Historically `codegen_function` accepted these as a `kwargs...` splat. Threading a
`kwargs...` bundle means every function that forwards codegen options gets specialized (and
recompiled) once per distinct *set* of keyword arguments, even when the values it forwards
are irrelevant to it. `CodegenFunctionOptions` is a single, non-parametric concrete type: no matter
which options are set, its type is always `CodegenFunctionOptions`, so functions that thread it
through only need to be compiled once.

The keyword constructor mirrors the historical keyword arguments of `codegen_function`
exactly, including their defaults. Unknown keyword arguments are ignored (as they were
previously silently dropped by `codegen_function`), preserving backwards compatibility.

Fields:

- `nanmath`: rewrite math functions to their `NaN`-safe variants.
- `wrap_code`: a 2-tuple of transformations applied to the out-of-place and in-place
  `Func`s respectively.
- `checkbounds`: emit bounds checks (when `false`, generated code is wrapped in
  `@inbounds`).
- `iip_config`: a 2-tuple of `Bool`s indicating whether to generate the out-of-place and
  in-place functions respectively.
- `sort_addmul`: sort the arguments of `+`/`*` for deterministic codegen.
- `optimize`: optimization rules to apply, or `nothing`.
- `similarto`: the array type (or name) to generate for array outputs, or `nothing` to
  infer it.
- `outputidxs`: explicit output indices for in-place array codegen, or `nothing`.
- `skipzeros`: skip writing structural zeros in in-place array codegen.
"""
struct CodegenFunctionOptions
    nanmath::Bool
    wrap_code::Tuple
    checkbounds::Bool
    iip_config::NTuple{2, Bool}
    sort_addmul::Bool
    optimize::Any
    similarto::Any
    outputidxs::Any
    skipzeros::Bool
end

function CodegenFunctionOptions(;
        nanmath::Bool = true, wrap_code::Tuple = (identity, identity),
        checkbounds = false, iip_config::NTuple{2, Bool} = (true, true),
        sort_addmul = false, optimize = nothing, similarto = nothing,
        outputidxs = nothing, skipzeros = false, kwargs...
    )
    return CodegenFunctionOptions(
        nanmath, wrap_code, checkbounds, iip_config, sort_addmul, optimize,
        similarto, outputidxs, skipzeros
    )
end

"""
    codegen_function(ir, expr, args[, options])
    codegen_function(ir, expr, args; kwargs...)

Generate out-of-place and in-place Julia function expressions for `expr` using a
`SymbolicUtils.IRStructure`. `args` is a vector whose elements describe each generated function
argument as either a scalar symbolic value or a collection of symbolic values.

Pass a [`CodegenFunctionOptions`](@ref) as `options`, or pass the same settings as keyword
arguments. The result is a 2-tuple containing the out-of-place and in-place function
expressions. The requested `iip_config` and the shape of `expr` determine which variants are
implemented.

This is the low-level code-generation interface for consumers that already maintain a
`SymbolicUtils.IRStructure`. Prefer [`build_function`](@ref) when starting from ordinary symbolic
expressions.
"""
function codegen_function(ir::IRStructure{VartypeT}, expr, args::Vector; kwargs...)
    return codegen_function(ir, expr, args, CodegenFunctionOptions(; kwargs...))
end

function codegen_function(
        ir::IRStructure{VartypeT}, expr, args::Vector, opts::CodegenFunctionOptions
    )
    (; nanmath, wrap_code, checkbounds, iip_config, sort_addmul, optimize) = opts
    args = canonicalize_args(args, !checkbounds)
    rewrites = Dict()
    if nanmath
        rewrites[:nanmath] = true
    end
    rewrites[:sort_addmul] = sort_addmul

    ir, expr = Code.apply_optimization_rules(ir, expr, optimize)

    if iip_config[1]
        oopfn = wrap_code[1](Func(args, [], expr))
        oop = Code.fast_toexpr(oopfn, ir, rewrites)
        if !checkbounds
            @assert Meta.isexpr(oop, :function)
            oop.args[2] = Expr(
                :macrocall, nameof(var"@inbounds"), LineNumberNode(0),
                Expr(:block, oop.args[2])
            )
        end
    else
        oop = get_unimplemented_fn(length(args), "out-of-place")
    end
    if iip_config[2] && SU.is_array_shape(SU.shape(expr))
        expr = SConst(expr)
        iipexpr = if Code.supports_with_allocator(expr)
            Code.with_allocator(IIP_ALLOCATOR, expr)
        else
            SU.Term{VartypeT}(
                copyto!, SArgsT((IIP_OUTSYM, expr));
                type = SU.symtype(expr), shape = SU.shape(expr)
            )
        end
        iipfn = wrap_code[2](Func([IIP_OUTSYM; args], [], iipexpr))
        iip = Code.fast_toexpr(iipfn, ir, rewrites)
        if !checkbounds
            @assert Meta.isexpr(iip, :function)
            iip.args[2] = Expr(
                :macrocall, nameof(var"@inbounds"), LineNumberNode(0),
                Expr(:block, iip.args[2])
            )
        end
    else
        iip = get_unimplemented_fn(length(args) + 1, "in-place")
    end
    return oop, iip
end

function codegen_function(
        ir::IRStructure{VartypeT}, expr::Union{Arr, Num, CallAndWrap, SymStruct},
        args::Vector, opts::CodegenFunctionOptions
    )
    return codegen_function(ir, unwrap(expr), args, opts)
end

function codegen_function(
        ir::IRStructure{VartypeT}, expr::AbstractArray, args::Vector,
        opts::CodegenFunctionOptions
    )
    (; similarto, nanmath, wrap_code, iip_config, outputidxs, skipzeros, checkbounds,
        optimize, sort_addmul) = opts
    args = canonicalize_args(args, !checkbounds)
    rewrites = Dict()
    if nanmath
        rewrites[:nanmath] = true
    end
    rewrites[:sort_addmul] = sort_addmul

    expr = _recursive_unwrap(expr)

    ir, expr = Code.apply_optimization_rules(ir, expr, optimize)

    i = findfirst(x -> x isa DestructuredArgs, args)
    if similarto === nothing
        similarto = i === nothing ? Array : (args[i]::DestructuredArgs).name
    end
    if iip_config[1]
        oopfn = wrap_code[1](Func(args, [], make_array(nothing, args, expr, similarto)))
        oop = Code.fast_toexpr(oopfn, ir, rewrites)
        if !checkbounds
            @assert Meta.isexpr(oop, :function)
            oop.args[2] = Expr(
                :macrocall, nameof(var"@inbounds"), LineNumberNode(0),
                Expr(:block, oop.args[2])
            )
        end
    else
        oop = get_unimplemented_fn(length(args), "out-of-place")
    end

    if iip_config[2]
        iipfn = wrap_code[2](
            Func(
                [IIP_OUTSYM; args], [], set_array(
                    nothing,
                    args,
                    IIP_OUTSYM,
                    outputidxs,
                    expr,
                    checkbounds,
                    skipzeros
                )
            )
        )
        iip = Code.fast_toexpr(iipfn, ir, rewrites)
        if !checkbounds
            @assert Meta.isexpr(iip, :function)
            iip.args[2] = Expr(
                :macrocall, nameof(var"@inbounds"), LineNumberNode(0),
                Expr(:block, iip.args[2])
            )
        end
    else
        iip = get_unimplemented_fn(length(args) + 1, "in-place")
    end
    return oop, iip
end
