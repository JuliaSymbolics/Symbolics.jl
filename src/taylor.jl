# TODO: around arbitrary point x = x0 ≠ 0
# TODO: error if x is not a "pure variable"
# TODO: optimize for multiple orders with loop/recursion
# TODO: get rational coefficients, not floats
"""
    taylor(f, x, n)

Calculate the `n`-th order term(s) in the Taylor series of the expression `f(x)` around `x = 0`.

Examples
========
```julia
julia> @variables x
1-element Vector{Num}:
 x

julia> taylor(exp(x), x, 0:3)
1.0 + x + 0.5(x^2) + 0.16666666666666666(x^3)
```
"""
function taylor(f, x, n::Int)
    D = Differential(x)
    n! = factorial(n)
    c = (D^n)(f) / n!
    c = expand_derivatives(c)
    c = substitute(c, x => 0)
    return c * x^n
end
function taylor(f, x, n::AbstractArray{Int})
    return sum(taylor.(f, x, n))
end
