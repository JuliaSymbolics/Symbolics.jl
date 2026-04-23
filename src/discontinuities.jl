"""
    rootfunction(f)

Given a function `f` with a discontinuity or discontinuous derivative, return the rootfinding
function of `f`. The rootfinding function `g` takes the same arguments as `f`, and is such
that `f` can be described as a piecewise function based on the sign of `g`, where each piece
is continuous and has a continuous derivative. The pieces are obtained using
`left_continuous_function(f)` and `right_continuous_function(f)`.

More formally,
```julia
f(args...) = if g(args...) < 0
    left_continuous_function(f)(args...)
else
    right_continuous_function(f)(args...)
end
```

For example, if `f` is `max(x, y)`, the root function is `(x, y) -> x - y` with
`left_continuous_function` as `(x, y) -> y` and `right_continuous_function` as
`(x, y) -> x`.

See also: [`left_continuous_function`](@ref), [`right_continuous_function`](@ref).
"""
function rootfunction end

"""
    left_continuous_function(f)

Given a function `f` with a discontinuity or discontinuous derivative, return a function
taking the same arguments as `f` which is continuous and has a continuous derivative
when `rootfinding_function(f)` is negative.

See also: [`rootfunction`](@ref).
"""
function left_continuous_function end

"""
    right_continuous_function(f)

Given a function `f` with a discontinuity or discontinuous derivative, return a function
taking the same arguments as `f` which is continuous and has a continuous derivative
when `rootfinding_function(f)` is positive.

See also: [`rootfunction`](@ref).
"""
function right_continuous_function end

"""
    @register_discontinuity f(arg1, arg2, ...) root_expr left_expr right_expr

Utility macro to register functions with discontinuities. The function `f` with
arguments `arg1, arg2, ...` has a `rootfunction` of `root_expr`, a
`left_continuous_function` of `left_expr` and `right_continuous_function` of
`right_expr`. `root_expr`, `left_expr` and `right_expr` are all expressions in terms
of `arg1, arg2, ...`.

For example, `max(x, y)` can be registered as `@register_discontinuity max(x, y) x - y y x`.

See also: [`rootfunction`](@ref)
"""
macro register_discontinuity(f, root, left, right)
    Meta.isexpr(f, :call) || error("Expected function call as first argument")
    args = f.args[2:end]
    fn = esc(f.args[1])
    rootname = gensym(:root)
    rootfn = :(function $rootname($(args...))
        $root
    end)
    leftname = gensym(:left)
    leftfn = :(function $leftname($(args...))
        $left
    end)
    rightname = gensym(:right)
    rightfn = :(function $rightname($(args...))
        $right
    end)
    return quote
        $rootfn
        (::$typeof($rootfunction))(::$typeof($fn)) = $rootname
        $leftfn
        (::$typeof($left_continuous_function))(::$typeof($fn)) = $leftname
        $rightfn
        (::$typeof($right_continuous_function))(::$typeof($fn)) = $rightname
    end
end

# a triangle function which is zero when x is a multiple of period
function _triangle(x, period)
    x /= 2period
    abs(x + 1 // 4 - floor(x + 3 // 4)) - 1 // 2
end

@register_discontinuity abs(x) x -x x
# just needs a rootfind to hit the discontinuity
@register_discontinuity mod(x, y) _triangle(x, y) mod(x, y) mod(x, y)
@register_discontinuity rem(x, y) _triangle(x, y) rem(x, y) rem(x, y)
@register_discontinuity div(x, y) _triangle(x, y) div(x, y) div(x, y)
@register_discontinuity max(x, y) x - y y x
@register_discontinuity min(x, y) x - y x y
@register_discontinuity NaNMath.max(x, y) x - y y x
@register_discontinuity NaNMath.min(x, y) x - y x y
@register_discontinuity <(x, y) x - y true false
@register_discontinuity <=(x, y) y - x false true
@register_discontinuity >(x, y) y - x true false
@register_discontinuity >=(x, y) x - y false true

"""
    majorization_function(f)

Given a function `f`, return a majorization function `m` for `f`. The function `m` should
have the signature `m(k, args...)` where `args...` are the same arguments as `f`. `k` is a
`Real` value which acts as an approximation factor. For higher `k`, the function `m` should
more closely approximate `f` over the domain. A majorization function is such that
`m(k, args...) >= f(args...)` for all `args...` in the domain.
"""
function majorization_function end

"""
    minorization_function(f)

Given a function `f`, return a minorization function `m` for `f`. The function `m` should
have the signature `m(k, args...)` where `args...` are the same arguments as `f`. `k` is a
`Real` value which acts as an approximation factor. For higher `k`, the function `m` should
more closely approximate `f` over the domain. A minorization function is such that
`m(k, args...) <= f(args...)` for all `args...` in the domain.
"""
function minorization_function end

"""
    approximation_function(f)

Given a function `f`, return an approximation function `appr` for `f`. The function `appr` should
have the signature `appr(k, args...)` where `args..` are the same arguments as `f`. `k` is a
`Real` value acting as an approximation factor. For higher `k`, the function `appr` should more
closely approximate `f` over the domain. The function `appr` offers no guarantees other than
infinite differentiability over the domain. At any point in the domain, it may evaluate to a
value greater or less than the value returned by `f` for the same point.
"""
function approximation_function end

function _logsumexp(m, a, b)
    a > b ? (a + log1p(exp(m * (b - a))) / m) : (b + log1p(exp(m * (a - b))) / m)
end

function _approx_min(m, a, b)
    -_logsumexp(m, -a, -b)
end

function _approx_abs(m, a)
    a + log1p(exp(-2m * a)) / m
end

function _approx_ge(m, a, b)
    (tanh(m * (a - b)) + 1) / 2
end

_approx_le(m, a, b) = _approx_ge(m, b, a)

function _sigmoid(x)
    if x > 0
        return one(x) / (one(x) + exp(-x))
    else
        tmp = exp(x)
        return tmp / (1 + tmp)
    end
end

function _sigder(x)
    tmp = _sigmoid(x)
    return tmp * (1 - tmp)
end

function _approx_eq(m, a, b)
    # `5m` to try and make `_approx_eq` similarly steep for the same `m` as other
    # approximators.
    4_sigder(5m * (a - b))
end

approximation_function(::typeof(max)) = _logsumexp
majorization_function(::typeof(max)) = _logsumexp
approximation_function(::typeof(min)) = _approx_min
minorization_function(::typeof(min)) = _approx_min
approximation_function(::typeof(NaNMath.max)) = _logsumexp
majorization_function(::typeof(NaNMath.max)) = _logsumexp
approximation_function(::typeof(NaNMath.min)) = _approx_min
minorization_function(::typeof(NaNMath.min)) = _approx_min
approximation_function(::typeof(abs)) = _approx_abs
majorization_function(::typeof(abs)) = _approx_abs

approximation_function(::typeof(>=)) = _approx_ge
approximation_function(::typeof(>)) = _approx_ge
approximation_function(::typeof(<=)) = _approx_le
approximation_function(::typeof(<)) = _approx_le

approximation_function(::typeof(==)) = _approx_eq
