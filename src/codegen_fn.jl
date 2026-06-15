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

function codegen_function(
        ir::IRStructure{VartypeT}, expr, args::Vector;
        nanmath::Bool = true, wrap_code::Tuple = (identity, identity),
        checkbounds = false, iip_config::NTuple{2, Bool} = (true, true), sort_addmul = false,
        optimize = nothing, kwargs...
    )
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
        args::Vector; kwargs...
    )
    return codegen_function(ir, unwrap(expr), args; kwargs...)
end

function codegen_function(
        ir::IRStructure{VartypeT}, expr::AbstractArray, args::Vector;
        similarto = nothing, nanmath::Bool = true, wrap_code::Tuple = (identity, identity),
        iip_config::NTuple{2, Bool} = (true, true), outputidxs = nothing,
        skipzeros = false, checkbounds = false, optimize = nothing, sort_addmul = false, kwargs...
    )
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
