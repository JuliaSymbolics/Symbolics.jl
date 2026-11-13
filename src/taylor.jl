
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

function _fraction_parts(ex)
    iscall(ex) || return (ex, 1)
    op = operation(ex)
    args = arguments(ex)
    if op === (/)
        n1, d1 = _fraction_parts(args[1])
        n2, d2 = _fraction_parts(args[2])
        return n1 * d2, d1 * n2
    elseif op === (*)
        parts = _fraction_parts.(args)
        return prod(first, parts), prod(last, parts)
    elseif op === (+)
        parts = _fraction_parts.(args)
        numerator = sum(parts[i][1] * prod(parts[j][2] for j in eachindex(parts) if j != i) for i in eachindex(parts))
        return numerator, prod(last, parts)
    elseif op === (^)
        exponent = value(args[2])
        if exponent isa Integer
            numerator, denominator = _fraction_parts(args[1])
            return exponent >= 0 ? (numerator^exponent, denominator^exponent) :
                (denominator^(-exponent), numerator^(-exponent))
        end
    end
    return ex, 1
end

_contains_nonfinite_constant(x::Number) = !isfinite(x)
_contains_nonfinite_constant(x::Num) = _contains_nonfinite_constant(unwrap(x))
function _contains_nonfinite_constant(x::BasicSymbolic{VartypeT})
    if iscall(x)
        return any(_contains_nonfinite_constant, arguments(x))
    elseif SymbolicUtils.isconst(x)
        return _contains_nonfinite_constant(unwrap_const(x))
    end
    return false
end
_contains_nonfinite_constant(x) = false

_iszero_coeff(c) = c isa Number ? iszero(c) : isequal(c, 0)

function _series_coeff(f, x, n; rationalize, kwargs...)
    numerator, denominator = value.(_fraction_parts(unwrap(f)))
    isequal(denominator, 1) && error("Cannot compute a finite Taylor coefficient because the expression is not a quotient with a vanishing denominator")
    isequal(value(simplify(denominator)), 0) && error("Cannot compute a Taylor coefficient with a zero denominator")

    # Search the vanishing order on raw derivatives: `taylor_coeff` would evaluate
    # `factorial(order)`, which overflows within the search bound (at 13 on 32-bit).
    D = Differential(x)
    max_order = n + 20
    order = 0
    derivative = expand_derivatives(denominator)
    while _iszero_coeff(value(substitute_in_deriv(derivative, x => 0; fold = Val(true), kwargs...)))
        order == max_order && error("Could not find a nonzero denominator Taylor coefficient through order $max_order at $x = 0")
        derivative = expand_derivatives(D(derivative))
        order += 1
    end

    numerator_coeffs = taylor_coeff(numerator, x, 0:(n + order); rationalize, kwargs...)
    denominator_coeffs = taylor_coeff(denominator, x, order:(n + order); rationalize, kwargs...)
    all(_iszero_coeff ∘ value, numerator_coeffs[1:order]) ||
        error("Cannot compute the Taylor coefficient of an expression with a pole at $x = 0")

    quotient_coeffs = Any[]
    for p in 0:n
        coefficient = numerator_coeffs[order + p + 1]
        for j in 0:(p - 1)
            coefficient -= quotient_coeffs[j + 1] * denominator_coeffs[p - j + 1]
        end
        push!(quotient_coeffs, coefficient / denominator_coeffs[1])
    end
    return value(quotient_coeffs[n + 1])
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

        isempty(n) && return similar(n, SymbolicT)

        # Differentiate once per order up to max(n), instead of
        # recomputing (D^k)(f) from scratch for each k.
        D = Differential(x)
        max_order = maximum(n)
        needed = Set{Int}(n)
        coeffs_by_order = Dict{Int, Any}()

        # expr holds the k-th derivative of f, updated incrementally below
        expr = expand_derivatives(f)
        for k in 0:max_order
            if k in needed
                c = value(substitute_in_deriv(expr, x => 0; fold = Val(true), kwargs...))
                coeffs_by_order[k] = _finalize_taylor_coeff(c, f, x, k; rationalize, kwargs...)
            end
            k < max_order && (expr = expand_derivatives(D(expr)))
        end

        return map(k -> coeffs_by_order[k], n)
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
    c = (D^n)(f)
    c = expand_derivatives(c)
    c = value(substitute_in_deriv(c, x => 0; fold = Val(true), kwargs...))
    return _finalize_taylor_coeff(c, f, x, n; rationalize, kwargs...)
end

function _finalize_taylor_coeff(c, f, x, n; rationalize, kwargs...)
    n! = factorial(n)
    if _contains_nonfinite_constant(c)
        c = n! * _series_coeff(f, x, n; rationalize, kwargs...)
    end
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

    if ns isa Number
        return taylor_coeff(f, x, ns; rationalize, kwargs...) * x^ns
    else
        # any iterable of orders (tuple, generator, ...) takes the incremental array path
        orders = ns isa AbstractArray ? ns : collect(ns)
        coeffs = taylor_coeff(f, x, orders; rationalize, kwargs...)
        return sum(c * x^n for (c, n) in zip(coeffs, orders))
    end
end
function taylor(f, x, x0, n; rationalize=true, kwargs...)
    # 1) substitute dummy x′ = x - x0
    name = Symbol(nameof(x), "′") # e.g. Symbol("x′")
    x′ = only(@variables $name)
    f = substitute_in_deriv(f, x => x′ + x0; kwargs...)

    # 2) expand f around x′ = 0
    s = try
        taylor(f, x′, n; rationalize, kwargs...)
    catch err
        if err isa ErrorException && occursin("pole at $x′ = 0", err.msg)
            throw(ErrorException(replace(err.msg, "pole at $x′ = 0" => "pole at $x = $x0")))
        end
        rethrow()
    end

    # 3) substitute back x = x′ + x0
    return substitute_in_deriv(s, x′ => x - x0; kwargs...)
end
