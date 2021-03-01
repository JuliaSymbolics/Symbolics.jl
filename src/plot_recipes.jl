# Plot recipe plotting an expression as a function
# expression may have 1 (ex1) or 2 (ex2) free variables
# plot of ex:𝑹 → 𝐑: `plot(ex1, a, b)`
# parametric plot: `plot(ex1, ex1′, a, b)`
# surface plot: `surface(xs, ys, ex2)`
# contourplot: `contour(xs, ys, ex2)`.
# not plot([ex1, ex1′], a, b); use plot(ex1, a, b); plot!(ex1′)
@recipe f(::Type{T}, v::T) where {T<:Num} = build_function(v, Symbolics.get_variables(v)...; expression=Val{false})
