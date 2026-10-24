"""
    inverse(f)

Given a single-input single-output function `f`, return its inverse `g`. This requires
that `f` is bijective. If `inverse` is defined for a function, `left_inverse` and
`right_inverse` should return `inverse(f)`. `inverse(g)` should also be defined to
return `f`.

See also: [`left_inverse`](@ref), [`right_inverse`](@ref), [`@register_inverse`](@ref).
"""
function inverse end

"""
    left_inverse(f)

Given a single-input single-output function `f`, return its left inverse `g`. This
requires that `f` is injective. If `left_inverse` is defined for a function,
`right_inverse` and `inverse` must not be defined and should error. `right_inverse(g)`
should also be defined to return `f`.

See also: [`inverse`](@ref), [`right_inverse`](@ref), [`@register_inverse`](@ref).
"""
function left_inverse end

"""
    right_inverse(f)

Given a single-input single-output function `f`, return its right inverse `g`. This
requires that `f` is surjective. If `right_inverse` is defined for a function,
`left_inverse` and `inverse` must not be defined and should error. `left_inverse(g)`
should also be defined to return `f`.

See also [`inverse`](@ref), [`left_inverse`](@ref), [`@register_inverse`](@ref).
"""
function right_inverse end

"""
    @register_inverse f g
    @register_inverse f g left
    @register_inverse f g right

Mark `f` and `g` as inverses of each other. By default, assume that `f` and `g` are
bijective. Also defines `left_inverse` and `right_inverse` to call `inverse`. If the
third argument is `left`, assume that `f` is injective and `g` is its left inverse. If
the third argument is `right`, assume that `f` is surjective and `g` is its right
inverse.
"""
macro register_inverse(f, g, dir::QuoteNode = :(:both))
    dir = dir.value
    f = esc(f)
    g = esc(g)
    if dir == :both
        quote
            (::$typeof($inverse))(::$typeof($f)) = $g
            (::$typeof($inverse))(::$typeof($g)) = $f
            (::$typeof($left_inverse))(::$typeof($f)) = $(inverse)($f)
            (::$typeof($right_inverse))(::$typeof($f)) = $(inverse)($f)
            (::$typeof($left_inverse))(::$typeof($g)) = $(inverse)($g)
            (::$typeof($right_inverse))(::$typeof($g)) = $(inverse)($g)
        end
    elseif dir == :left
        quote
            (::$typeof($left_inverse))(::$typeof($f)) = $g
            (::$typeof($right_inverse))(::$typeof($g)) = $f
        end
    elseif dir == :right
        quote
            (::$typeof($right_inverse))(::$typeof($f)) = $g
            (::$typeof($left_inverse))(::$typeof($g)) = $f
        end
    else
        throw(ArgumentError("The third argument to `@register_inverse` must be `left` or `right`"))
    end
end

"""
    $(TYPEDSIGNATURES)

Check if the provided function has an inverse defined via [`inverse`](@ref). Uses
`hasmethod` to perform the check.
"""
has_inverse(::T) where {T} = hasmethod(inverse, Tuple{T})

"""
    $(TYPEDSIGNATURES)

Check if the provided function has a left inverse defined via [`left_inverse`](@ref)
Uses `hasmethod` to perform the check.
"""
has_left_inverse(::T) where {T} = hasmethod(left_inverse, Tuple{T})

"""
    $(TYPEDSIGNATURES)

Check if the provided function has a left inverse defined via [`left_inverse`](@ref)
Uses `hasmethod` to perform the check.
"""
has_right_inverse(::T) where {T} = hasmethod(right_inverse, Tuple{T})

"""
    $(TYPEDSIGNATURES)

A simple utility function which returns the square of the input. Used to define
the inverse of `sqrt`.
"""
square(x) = x ^ 2

"""
    $(TYPEDSIGNATURES)

A simple utility function which returns the cube of the input. Used to define
the inverse of `cbrt`.
"""
cube(x) = x ^ 3

"""
    $(TYPEDSIGNATURES)

A simple utility function which takes `x` and returns `acos(x) / pi`. Used to
define the inverse of `acospi`.
"""
acosbypi(x) = acos(x) / pi

@register_inverse sin asin
@register_inverse cos acos
@register_inverse tan atan
@register_inverse csc acsc
@register_inverse sec asec
@register_inverse cot acot
@register_inverse sind asind
@register_inverse cosd acosd
@register_inverse tand atand
@register_inverse cscd acscd
@register_inverse secd asecd
@register_inverse cotd acotd
@register_inverse sinh asinh
@register_inverse cosh acosh
@register_inverse tanh atanh
@register_inverse csch acsch
@register_inverse sech asech
@register_inverse coth acoth
@register_inverse cospi acosbypi
@register_inverse SpecialFunctions.digamma SpecialFunctions.invdigamma
@register_inverse log exp
@register_inverse log2 exp2
@register_inverse log10 exp10
@register_inverse log1p expm1
@register_inverse deg2rad rad2deg
@register_inverse sqrt square :left
@register_inverse cbrt cube
@register_inverse NaNMath.sin NaNMath.asin
@register_inverse NaNMath.cos NaNMath.acos
# can't use macro since it would be a re-definition of `inverse(atan)`
inverse(::typeof(NaNMath.tan)) = inverse(tan)
inverse(::typeof(NaNMath.acosh)) = inverse(acosh)
inverse(::typeof(NaNMath.atanh)) = inverse(atanh)
inverse(::typeof(NaNMath.log)) = inverse(log)
inverse(::typeof(NaNMath.log10)) = inverse(log10)
inverse(::typeof(NaNMath.log1p)) = inverse(log1p)
inverse(::typeof(NaNMath.log2)) = inverse(log2)
left_inverse(::typeof(NaNMath.sqrt)) = left_inverse(sqrt)
# inverses of solve helpers
left_inverse(::typeof(ssqrt)) = left_inverse(sqrt)
left_inverse(::typeof(scbrt)) = left_inverse(cbrt)
left_inverse(::typeof(slog)) = left_inverse(log)

function inverse(f::ComposedFunction)
    return inverse(f.inner) ∘ inverse(f.outer)
end
has_inverse(f::ComposedFunction) = has_inverse(f.inner) && has_inverse(f.outer)
function left_inverse(f::ComposedFunction)
    return left_inverse(f.inner) ∘ left_inverse(f.outer)
end
function has_left_inverse(f::ComposedFunction)
    return has_left_inverse(f.inner) && has_left_inverse(f.outer)
end
function right_inverse(f::ComposedFunction)
    return right_inverse(f.inner) ∘ right_inverse(f.outer)
end
function has_right_inverse(f::ComposedFunction)
    return has_right_inverse(f.inner) && has_right_inverse(f.outer)
end
