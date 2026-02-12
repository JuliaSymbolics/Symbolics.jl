module SymbolicsReactantExt

using Symbolics
using Symbolics: Arr, Num, JuliaTarget, ReactantTarget, SerialForm
using Reactant: Reactant, @code_hlo, @compile, @allowscalar

function __concrete_input(::Type{T}, x::Union{<:Array{Num}, Arr}) where {T}
    return Reactant.to_rarray(zeros(T, size(x)...))
end
function __concrete_input(::Type{T}, ::Num) where {T}
    return Reactant.to_rarray(zero(T); track_numbers = true)
end

function Symbolics._build_function(
        target::ReactantTarget, rhss::AbstractArray, args...;
        eltypes = nothing,
        expression = Val{true}, compile_options = Reactant.CompileOptions(), kwargs...
)
    _, fn_iip_ = build_function(
        rhss, args...; target = JuliaTarget(), expression = Val(false),
        kwargs..., parallel = SerialForm()
    )

    if eltypes === nothing
        eltypes = [Float32 for _ in 1:(length(args) + 1)]
    elseif eltypes isa Type
        eltypes = [eltypes for _ in 1:(length(args) + 1)]
    else
        @assert length(eltypes) == length(args) + 1
    end

    out_concrete = __concrete_input(eltypes[1], rhss)
    args_concrete = map(__concrete_input, eltypes[2:end], args)

    fn_iip = (du, u, p, t) -> begin
        @allowscalar fn_iip_(du, u, p, t)
    end

    fn_oop = (u, p, t) -> begin
        out = similar(u, size(rhss)...)
        fn_iip(out, u, p, t)
        return out
    end

    if expression == Val{true}
        oop_mlir = @code_hlo compile_options=compile_options fn_oop(args_concrete...)
        iip_mlir = @code_hlo compile_options=compile_options fn_iip(
            out_concrete, args_concrete...)
        return oop_mlir, iip_mlir
    else
        fn_oop_compiled = @compile compile_options=compile_options fn_oop(args_concrete...)
        fn_iip_compiled = @compile compile_options=compile_options fn_iip(
            out_concrete, args_concrete...)
        return fn_oop_compiled, fn_iip_compiled
    end
end

end
