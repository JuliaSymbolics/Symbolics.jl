"""
    series(cs, x, [x0=0,], ns=0:length(cs)-1)

Return the power series in `x` around `x0` to the powers `ns` with coefficients `cs`.

    series(y, x, [x0=0,] ns)

Return the power series in `x` around `x0` to the powers `ns` with coefficients automatically created from the variable `y`.

Examples
========

```julia
julia> @variables x y[0:3] z
3-element Vector{Any}:
 x
  y[0:3]
 z

julia> series(y, x, 2)
y[0] + (-2 + x)*y[1] + ((-2 + x)^2)*y[2] + ((-2 + x)^3)*y[3]

julia> series(z, x, 2, 0:3)
z[0] + (-2 + x)*z[1] + ((-2 + x)^2)*z[2] + ((-2 + x)^3)*z[3]
```
"""
function series(cs::AbstractArray, x::Number, x0::Number, ns::AbstractArray = 0:length(cs)-1)
    length(cs) == length(ns) || error("There are different numbers of coefficients and orders")
    s = sum(c * (x - x0)^n for (c, n) in zip(cs, ns))
    return s
end
function series(cs::AbstractArray, x::Number, ns::AbstractArray = 0:length(cs)-1)
    return series(cs, x, 0, ns)
end
function series(y::Num, x::Number, x0::Number, ns::AbstractArray)
    cs, = @variables $(nameof(y))[ns]
    return series(cs, x, x0, ns)
end
function series(y::Num, x::Number, ns::AbstractArray)
    return series(y, x, 0, ns)
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
        c = Base.rationalize(c) # convert integers/floats to rational numbers; avoid name clash between rationalize and Base.rationalize()
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
