
# Alex: make sure `Num`s are not processed here as they'd break it.
_is_const_number(x::Number) = true
function _is_const_number(x::SymbolicUtils.BasicSymbolic)
    !iscall(x) && return false
    all(_is_const_number, arguments(x))
end

_postprocess_root(x) = x

function _postprocess_root(x::Number)
    # N // 1 => N
    if x isa Rational
        if isone(denominator(x))
            return numerator(x)
        end
    end

    # N//1 + M//1*im => N + M*im
    if x isa Complex{Rational{T}} where {T <: Integer}
        if isone(denominator(real(x))) && isone(denominator(imag(x)))
            return numerator(real(x)) + numerator(imag(x))*im
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

    x = Symbolics.term(operation(x), map(_postprocess_root, arguments(x))...)

    # sqrt(0), cbrt(0) => 0
    # sqrt(1), cbrt(1) => 1
    if iscall(x) &&
       (operation(x) === sqrt || operation(x) === cbrt || operation(x) === ssqrt ||
        operation(x) === scbrt)
        arg = arguments(x)[1]
        if isequal(arg, 0) || isequal(arg, 1)
            return arg
        end
    end

    # (X)^0 => 1
    if iscall(x) && operation(x) === (^) && isequal(arguments(x)[2], 0)
        return 1
    end

    # (X)^1 => X
    if iscall(x) && operation(x) === (^) && isequal(arguments(x)[2], 1)
        return arguments(x)[1]
    end

    # sqrt((N / D)^2 * M) => N / D * sqrt(M)
    if iscall(x) && (operation(x) === sqrt || operation(x) === ssqrt)
        function squarefree_decomp(x::Integer)
            square, squarefree = big(1), big(1)
            for (p, d) in collect(Primes.factor(abs(x)))
                q, r = divrem(d, 2)
                square *= p^q
                squarefree *= p^r
            end
            square, squarefree
        end
        arg = arguments(x)[1]
        if arg isa Integer
            square, squarefree = squarefree_decomp(arg)
            if arg < 0
                square = im * square
            end
            if !isone(square)
                return square * Symbolics.term(Symbolics.operation(x), squarefree)
            end
        elseif arg isa Rational
            n, d = numerator(arg), denominator(arg)
            n_square, n_squarefree = squarefree_decomp(n)
            if n < 0
                n_square = im * n_square
            end
            d_square, d_squarefree = squarefree_decomp(d)
            nd_square = n_square // d_square
            nd_squarefree = n_squarefree // d_squarefree
            if !isone(nd_square)
                return nd_square * Symbolics.term(Symbolics.operation(x), nd_squarefree)
            end
        end
    end

    # (sqrt(N))^M => N^div(M, 2)*sqrt(N)^(mod(M, 2))
    if iscall(x) && operation(x) === (^)
        arg1, arg2 = arguments(x)
        if iscall(arg1) && (operation(arg1) === sqrt || operation(arg1) === ssqrt)
            if arg2 isa Integer
                isequal(arg2, 2) && return arguments(arg1)[1]
                q, r = divrem(arg2, 2)
                if isequal(r, 1)
                    return arguments(arg1)[1]^q * arg1
                else
                    return arguments(arg1)[1]^q
                end
            end
        end
    end

    return x
end

function postprocess_root(x)
    math_consts = (Base.MathConstants.pi, Base.MathConstants.e)
    while true
        old_x = deepcopy(x)
        contains_math_const = any([Symbolics.n_occurrences(x, c) > 0 for c in math_consts])
        if contains_math_const
            x = _postprocess_root(x)
        else
            x = x |> expand |> _postprocess_root |> expand
        end
        isequal(typeof(old_x), typeof(x)) && isequal(old_x, x) && return x
    end
    x # unreachable
end
