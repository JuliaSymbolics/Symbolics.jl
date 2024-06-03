"""
    SymbolicsSparsityDetector <: ADTypes.AbstractSparsityDetector

Sparsity detection algorithm based on the [Symbolics.jl tracing system](https://symbolics.juliasymbolics.org/stable/manual/sparsity_detection/).

This type makes Symbolics.jl compatible with the [ADTypes.jl sparsity detection framework](https://sciml.github.io/ADTypes.jl/stable/#Sparsity-detector). The following functions are implemented:

- `ADTypes.jacobian_sparsity` based on [`Symbolics.jacobian_sparsity`](@ref)
- `ADTypes.hessian_sparsity` based on [`Symbolics.hessian_sparsity`](@ref)

# Reference

> [Sparsity Programming: Automated Sparsity-Aware Optimizations in Differentiable Programming](https://openreview.net/forum?id=rJlPdcY38B), Gowda et al. (2019)
"""
struct SymbolicsSparsityDetector <: ADTypes.AbstractSparsityDetector end

function ADTypes.jacobian_sparsity(f, x::AbstractArray, ::SymbolicsSparsityDetector)
    y = similar(f(x))
    f!(y, x) = copyto!(y, f(x))
    return jacobian_sparsity(f!, y, x)
end

function ADTypes.jacobian_sparsity(f!, y::AbstractArray, x::AbstractArray, ::SymbolicsSparsityDetector)
    f!_vec(y_vec, x_vec) = f!(reshape(y_vec, size(y)), reshape(x_vec, size(x)))
    return jacobian_sparsity(f!_vec, vec(y), vec(x))
end

function ADTypes.hessian_sparsity(f, x::AbstractArray, ::SymbolicsSparsityDetector)
    f_vec(x_vec) = f(reshape(x_vec, size(x)))
    return hessian_sparsity(f_vec, vec(x))
end
