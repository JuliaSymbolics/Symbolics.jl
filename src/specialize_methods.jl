"""
    specialize_methods(func, abstract_arg_types, inner_func, mods=nothing)

For any method that implements `func` with signature 
fitting `abstract_arg_types`, define methods for corresponding 
symbolic types that pass all arguments to `inner_func`.
`mods` is an optional list of modules to look for methods in.
"""
function specialize_methods(func, abstract_arg_types, inner_func, mods=nothing)
    ms = isnothing(mods) ? methods(func, abstract_arg_types) : methods(func, abstract_arg_types, mods)
    for m in ms
        mod = m.module
        if mod != @__MODULE__ # do not overwrite method definitions from within this module itself, else: precompilation warnings
            sig = ExprTools.signature(m; extra_hygiene=true)
            fname = sig[:name]
            args = sig[:args]
            kwargs = get(sig, :kwargs, Symbol[])
            whereparams = get(sig, :whereparams, Symbol[])
            args_names = expr_argname.(args)
            kwargs_names = expr_kwargname.(kwargs)
            body = :($(inner_func)($(args_names...); $(kwargs_names...)))
            Base.eval(
                @__MODULE__, 
                wrap_func_expr(
                    mod, fname, args, kwargs, args_names, kwargs_names, whereparams, body;
                    abstract_arg_types
                )
            )
        end#of `mod != @__MODULE__`
    end#of `for m in ms`
end

function specialize_methods(mods=nothing)
    specialize_methods(Base.:(*), (AbstractMatrix, AbstractVector), _matvec, mods)
    specialize_methods(Base.:(*), (AbstractMatrix, AbstractMatrix), _matmul, mods)
end