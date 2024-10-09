"""
    series(y, x, [x0=0,] ns; name = nameof(y))

Expand the variable `y` in a power series in the variable `x` around `x0` to orders `ns`.

Examples
========

```julia
julia> @variables z ϵ
2-element Vector{Num}:
 z
 ϵ

julia> series(z, ϵ, 2, 0:3)
z[0] + z[1]*(-2 + ϵ) + z[2]*((-2 + ϵ)^2) + z[3]*((-2 + ϵ)^3)
```
"""
function series(y, x, x0, ns; name = nameof(y))
    c, = @variables $name[ns]
    return sum(c[n] * (x - x0)^n for n in ns)
end
function series(y, x, ns; kwargs...)
    return series(y, x, 0, ns; kwargs...)
end

"""
    taylor_coeff(f, x[, n]; rationalize=true)

Calculate the `n`-th order coefficient(s) in the Taylor series of `f` around `x = 0`.

Examples
========
```julia
julia> @variables x y
2-element Vector{Num}:
 x
 y

julia> taylor_coeff(series(y, x, 0:5), x, 0:2:4)
3-element Vector{Num}:
 y[0]
 y[2]
 y[4]
```
"""
function taylor_coeff(f, x, n = missing; rationalize=true)
    if n isa AbstractArray
        # return array of expressions/equations for each order
        return taylor_coeff.(Ref(f), Ref(x), n; rationalize)
    elseif f isa Equation
        if ismissing(n)
            # assume user wants maximum order in the equation
            n = 0:max(degree(f.lhs, x), degree(f.rhs, x))
            return taylor_coeff(f, x, n; rationalize)
        else
            # return new equation with coefficients of each side
            return taylor_coeff(f.lhs, x, n; rationalize) ~ taylor_coeff(f.rhs, x, n; rationalize)
        end
    elseif ismissing(n)
        # assume user wants maximum order in the expression
        n = 0:degree(f, x)
        return taylor_coeff(f, x, n; rationalize)
    end

    # TODO: error if x is not a "pure variable"
    D = Differential(x)
    n! = factorial(n)
    c = (D^n)(f) / n! # TODO: optimize the implementation for multiple n with a loop that avoids re-differentiating the same expressions
    c = expand_derivatives(c)
    c = substitute(c, x => 0)
    if rationalize && unwrap(c) isa Number
        # TODO: make rational coefficients "organically" and not using rationalize (see https://github.com/JuliaSymbolics/Symbolics.jl/issues/1299)
        c = unwrap(c)
        c = isinteger(c) ? Int(c) : Base.rationalize(c) # convert non-integer floats to rational numbers; avoid name clash between rationalize and Base.rationalize()
    end
    return c
end

"""
    taylor(f, x, [x0=0,] n; rationalize=true)

Calculate the `n`-th order term(s) in the Taylor series of `f` around `x = x0`.
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

julia> taylor(√(x), x, 1, 0:3)
1 + (1//2)*(-1 + x) - (1//8)*((-1 + x)^2) + (1//16)*((-1 + x)^3)

julia> isequal(taylor(exp(im*x), x, 0:5), taylor(exp(im*x), x, 0:5))
true
```
"""
function taylor(f, x, ns; kwargs...)
    if f isa AbstractArray
        return taylor.(f, Ref(x), Ref(ns); kwargs...)
    elseif f isa Equation
        return taylor(f.lhs, x, ns; kwargs...) ~ taylor(f.rhs, x, ns; kwargs...)
    end

    return sum(taylor_coeff(f, x, n; kwargs...) * x^n for n in ns)
end
function taylor(f, x, x0, n; kwargs...)
    # 1) substitute dummy x′ = x - x0
    name = Symbol(nameof(x), "′") # e.g. Symbol("x′")
    x′ = only(@variables $name)
    f = substitute(f, x => x′ + x0)

    # 2) expand f around x′ = 0
    s = taylor(f, x′, n; kwargs...)

    # 3) substitute back x = x′ + x0
    return substitute(s, x′ => x - x0)
end
