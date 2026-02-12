module SymbolicsReactantExt

using Symbolics
using Symbolics: JuliaTarget, ReactantTarget, SerialForm
using Reactant: Reactant, @code_hlo, @compile, @allowscalar

function Symbolics._build_function(target::ReactantTarget, rhss::AbstractArray, args...;
        expression = Val{true}, kwargs...)
    _, fn_iip_ = build_function(
        rhss, args...; target = JuliaTarget(), expression = Val(false),
        kwargs..., parallel = SerialForm()
    )

    # TODO: allow configuring the eltypes
    out_concrete = Reactant.to_rarray(zeros(Float32, size(rhss)...))
    u_concrete = Reactant.to_rarray(zeros(Float32, size(args[1])...))
    p_concrete = Reactant.to_rarray(zeros(Float32, size(args[2])...))
    t_concrete = Reactant.to_rarray(zero(Float32); track_numbers = true)

    fn_iip = (du, u, p, t) -> begin
        @allowscalar fn_iip_(du, u, p, t)
    end

    fn_oop = (u, p, t) -> begin
        out = similar(u, size(rhss)...)
        fn_iip(out, u, p, t)
        return out
    end

    compile_options = target.compile_options === nothing ? Reactant.CompileOptions() :
                      target.compile_options

    if expression == Val{true}
        oop_mlir = @code_hlo compile_options=compile_options fn_oop(
            u_concrete, p_concrete, t_concrete)
        iip_mlir = @code_hlo compile_options=compile_options fn_iip(
            out_concrete, u_concrete, p_concrete, t_concrete)
        return oop_mlir, iip_mlir
    else
        fn_oop_compiled = @compile compile_options=compile_options fn_oop(
            u_concrete, p_concrete, t_concrete)
        fn_iip_compiled = @compile compile_options=compile_options fn_iip(
            out_concrete, u_concrete, p_concrete, t_concrete)
        return fn_oop_compiled, fn_iip_compiled
    end
end

end
