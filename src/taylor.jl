# TODO: around arbitrary point x = x0 ≠ 0
# TODO: error if x is not a "pure variable"
# TODO: optimize for multiple orders with loop/recursion
"""
    taylor(f, x, n; rationalize=true)

Calculate the `n`-th order term(s) in the Taylor series of the expression `f(x)` around `x = 0`.
If `rationalize`, float coefficients are approximated as rational numbers (this can produce unexpected results for irrational numbers, for example).

Examples
========
```julia
julia> @variables x
1-element Vector{Num}:
 x

julia> taylor(exp(x), x, 0:3)
1 + x + (1//2)*(x^2) + (1//6)*(x^3)

julia> taylor(exp(x), x, 0:3; rationalize=false)
1.0 + x + 0.5(x^2) + 0.16666666666666666(x^3)
```
"""
function taylor(f, x, n::Int; rationalize=true)
    D = Differential(x)
    n! = factorial(n)
    c = (D^n)(f) / n!
    c = expand_derivatives(c)
    c = substitute(c, x => 0)
    if rationalize && unwrap(c) isa Number
        # TODO: make rational coefficients "organically" and not using rationalize (see https://github.com/JuliaSymbolics/Symbolics.jl/issues/1299)
        c = unwrap(c)
        c = c % 1.0 == 0.0 ? Int(c) : Base.rationalize(c) # avoid nameclash with keyword argument
    end
    return c * x^n
end
function taylor(f, x, n::AbstractArray{Int}; kwargs...)
    return sum(taylor.(f, x, n; kwargs...))
end
