
"""
    series(cs, x, [x0=0,], ns=0:length(cs)-1)

Return the power series in `x` around `x0` to the powers `ns` with coefficients `cs`.

    series(y, x, [x0=0,] ns)

Return the power series in `x` around `x0` to the powers `ns` with coefficients automatically created from the variable `y`.

Examples
========

```jldoctest
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
    taylor_coeff(f, x[, n]; rationalize=true, kwargs...)

Calculate the `n`-th order coefficient(s) in the Taylor series of `f` around `x = 0`.
If `rationalize`, float coefficients are approximated as rational numbers (this can produce unexpected results for irrational numbers, for example).
Keyword arguments `kwargs...` are forwarded to internal `substitute()` calls.

Examples
========
```jldoctest
julia> @variables x y
2-element Vector{Num}:
 x
 y

julia> taylor_coeff(series(y, x, 0:5), x, 0:2:4)
3-element Vector{SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymReal}}:
 y[0]
 y[2]
 y[4]
```
"""
function taylor_coeff(f, x, n = missing; rationalize=true, kwargs...)
    if n isa AbstractArray
        if f isa Equation
            # apply incremental computation to each side of the equation
            lhs_coeffs = taylor_coeff(f.lhs, x, n; rationalize, kwargs...)
            rhs_coeffs = taylor_coeff(f.rhs, x, n; rationalize, kwargs...)
            return Equation[l ~ r for (l, r) in zip(lhs_coeffs, rhs_coeffs)]
        end

        # Differentiate once per order up to max(n), instead of
        # recomputing (D^k)(f) from scratch for each k.
        coeffs = similar(n, SymbolicT)
        isempty(n) && return coeffs

        D = Differential(x)
        max_order = maximum(n)
        IndexT = eltype(eachindex(n))
        indices = Dict{Int, Vector{IndexT}}()
        for i in eachindex(n)
            push!(get!(Vector{IndexT}, indices, n[i]), i)
        end

        # expr holds the k-th derivative of f, updated incrementally below
        expr = expand_derivatives(f)
        for k in 0:max_order
            if haskey(indices, k)
                c = value(substitute_in_deriv(expr, x => 0; fold = Val(true), kwargs...))
                k! = factorial(k)
                if !(c isa BasicSymbolic{VartypeT}) && isinteger(c)
                    c = Integer(c)
                    c //= k!
                elseif c isa Rational
                    c //= k!
                else
                    c /= k!
                end
                if rationalize && unwrap(c) isa Number
                    c = unwrap(c)
                    c = Base.rationalize(c)
                end
                if !(c isa BasicSymbolic{VartypeT})
                    c = SConst(c)
                end
                for i in indices[k]
                    coeffs[i] = c
                end
            end
            k < max_order && (expr = expand_derivatives(D(expr)))
        end

        return coeffs
    elseif f isa Equation
        if ismissing(n)
            # assume user wants maximum order in the equation
            n = 0:max(degree(f.lhs, x), degree(f.rhs, x))
            return taylor_coeff(f, x, n; rationalize, kwargs...)
        else
            # return new equation with coefficients of each side
            return taylor_coeff(f.lhs, x, n; rationalize, kwargs...) ~ taylor_coeff(f.rhs, x, n; rationalize, kwargs...)
        end
    elseif ismissing(n)
        # assume user wants maximum order in the expression
        n = 0:degree(f, x)
        return taylor_coeff(f, x, n; rationalize, kwargs...)
    end

    # TODO: error if x is not a "pure variable"
    D = Differential(x)
    n! = factorial(n)
    c = (D^n)(f)
    c = expand_derivatives(c)
    c = value(substitute_in_deriv(c, x => 0; fold = Val(true), kwargs...))
    if !(c isa BasicSymbolic{VartypeT}) && isinteger(c)
        c = Integer(c)
        c //= n!
    elseif c isa Rational
        c //= n!
    else
        c /= n!
    end
    if rationalize && unwrap(c) isa Number
        # TODO: make rational coefficients "organically" and not using rationalize (see https://github.com/JuliaSymbolics/Symbolics.jl/issues/1299)
        c = unwrap(c)
        c = Base.rationalize(c) # convert integers/floats to rational numbers; avoid name clash between rationalize and Base.rationalize()
    end
    return c
end

"""
    taylor(f, x, [x0=0,] n; rationalize=true, kwargs...)

Calculate the `n`-th order term(s) in the Taylor series of `f` around `x = x0`.
If `rationalize`, float coefficients are approximated as rational numbers (this can produce unexpected results for irrational numbers, for example).
Keyword arguments `kwargs...` are forwarded to internal `substitute()` calls.

Examples
========
```jldoctest
julia> @variables x
1-element Vector{Num}:
 x

julia> taylor(exp(x), x, 0:3)
1 + x + (1//2)*(x^2) + (1//6)*(x^3)

julia> taylor(exp(x), x, 0:3; rationalize=false)
1 + x + (1//2)*(x^2) + (1//6)*(x^3)

julia> taylor(√(x), x, 1, 0:3)
1 + (1//2)*(-1 + x) - (1//8)*((-1 + x)^2) + (1//16)*((-1 + x)^3)

julia> isequal(taylor(exp(im*x), x, 0:5), taylor(exp(im*x), x, 0:5))
true
```
"""
function taylor(f, x, ns; rationalize=true, kwargs...)
    if f isa AbstractArray
        return taylor.(f, Ref(x), Ref(ns); rationalize, kwargs...)
    elseif f isa Equation
        return taylor(f.lhs, x, ns; rationalize, kwargs...) ~ taylor(f.rhs, x, ns; rationalize, kwargs...)
    end

    if ns isa AbstractArray
        # Use the incremental path: compute all coefficients at once, then sum
        coeffs = taylor_coeff(f, x, ns; rationalize, kwargs...)
        return sum(c * x^n for (c, n) in zip(coeffs, ns))
    else
        return taylor_coeff(f, x, ns; rationalize, kwargs...) * x^ns
    end
end
function taylor(f, x, x0, n; rationalize=true, kwargs...)
    # 1) substitute dummy x′ = x - x0
    name = Symbol(nameof(x), "′") # e.g. Symbol("x′")
    x′ = only(@variables $name)
    f = substitute_in_deriv(f, x => x′ + x0; kwargs...)

    # 2) expand f around x′ = 0
    s = taylor(f, x′, n; rationalize, kwargs...)

    # 3) substitute back x = x′ + x0
    return substitute_in_deriv(s, x′ => x - x0; kwargs...)
end
