# import SymbolicUtils
# import 
import TermInterface: maketerm
# 
# using Test

# Alex: can be used in the following way.
#
# @variables x
# sol = RootFinding.solve(x^12 - 1, x)
# map(RootFinding.postprocess_root, sol)

# Alex: make sure `Num`s are not processed here as they'd break it.
_is_const_number(x::Number) = true
function _is_const_number(x::SymbolicUtils.BasicSymbolic)
    !iscall(x) && return false
    all(_is_const_number, arguments(x))
end

SymbolicUtils.@syms __x
# @test !_is_const_number(__x) && !_is_const_number(sqrt(__x))
# @test _is_const_number(1) && _is_const_number(2 // 3) && _is_const_number(3 + 4im)
# @test _is_const_number(SymbolicUtils.term(sqrt, 2) + 21)
# @test _is_const_number((SymbolicUtils.term(exp, 2) * SymbolicUtils.term(exp, 2)) // 99)

_postprocess_root(x) = x

function _postprocess_root(x::Number)
    # N // 1 => N
    if x isa Rational
        if isone(denominator(x))
            return big(numerator(x))
        end
    end

    # A + im*0 => A
    if x isa Complex
        if iszero(imag(x))
            # Alex: maybe return _postprocess_root(real(x))
            return real(x)
        end
    end

    return x
end

function _postprocess_root(x::SymbolicUtils.BasicSymbolic)
    !iscall(x) && return x

    x = maketerm(
        typeof(x),
        operation(x),
        map(_postprocess_root, arguments(x)),
        nothing
    )

    # sqrt(0), cbrt(0) => 0
    # sqrt(1), cbrt(1) => 1
    if iscall(x) && (operation(x) === sqrt || operation(x) === cbrt || operation(x) === ssqrt || operation(x) === scbrt)
        arg = arguments(x)[1]
        if isequal(arg, 0) || isequal(arg, 1)
            return arg
        end
    end

    # sqrt(N^2) => N
    if iscall(x) && (operation(x) === sqrt || operation(x) === ssqrt) 
        arg = arguments(x)[1]
        if arg isa Integer && arg > 0 && arg == (isqrt(arg))^2
            return isqrt(arg)
        elseif arg isa Integer && arg < 0 && -arg == (isqrt(-arg))^2
            return im*isqrt(-arg)
        end
    end

    
    # (sqrt(N))^2 => N
    if iscall(x) && operation(x) === (^) && isequal(arguments(x)[2], 2)
        arg1 = arguments(x)[1]
        if iscall(arg1) && (operation(arg1) === sqrt || operation(arg1) === ssqrt)
            return arguments(arg1)[1]
        end
    end

    return x
end

function postprocess_root(x)
    _postprocess_root(x)
end



# some tests
# SymbolicUtils.@syms __x
# __symsqrt(x) = SymbolicUtils.term(ssqrt, x)
# @test postprocess_root(2 // 1) == 2 && postprocess_root(2 + 0*im) == 2
# @test postprocess_root(__symsqrt(__symsqrt(0)) - 11) == -11
# @test postprocess_root(3*__symsqrt(2)^2) == 6
# @test postprocess_root(__symsqrt(4)) == 2
# @test isequal(postprocess_root(__symsqrt(__x)^2), __symsqrt(__x)^2)

