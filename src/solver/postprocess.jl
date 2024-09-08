# Alex: make sure `Num`s are not processed here as they'd break it.
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
    oper = operation(x)

    # sqrt(0), cbrt(0) => 0
    # sqrt(1), cbrt(1) => 1
    if (oper === sqrt || oper === cbrt || oper === ssqrt ||
        oper === scbrt)
        arg = arguments(x)[1]
        if isequal(arg, 0) || isequal(arg, 1)
            return arg
        end
    end

    # (X)^0 => 1
    if oper === (^) && isequal(arguments(x)[2], 0)
        return 1
    end

    # (X)^1 => X
    if oper === (^) && isequal(arguments(x)[2], 1)
        return arguments(x)[1]
    end

    # sqrt((N / D)^2 * M) => N / D * sqrt(M)
    if (oper === sqrt || oper === ssqrt)
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
    if oper === (^)
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

    trig_simplified = check_trig_consts(x)
    !isequal(trig_simplified, x) && return trig_simplified

    if oper === (+)
        args = arguments(x)
        for arg in args
            if isequal(arg, 0)
                after_removing = setdiff(args, arg)
                isone(length(after_removing)) && return after_removing[1]
                return Symbolics.term(+, after_removing)
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

function check_trig_consts(x)
    !iscall(x) && return x

    oper = operation(x)
    inv_opers = [asin, acos, atan]
    inv_exacts = [0, Symbolics.term(*, pi),
        Symbolics.term(/,pi,3),
        Symbolics.term(/, pi, 2),
        Symbolics.term(/, Symbolics.term(*, 2, pi), 3),
        Symbolics.term(/, pi, 6),
        Symbolics.term(/, Symbolics.term(*, 5, pi), 6),
        Symbolics.term(/, pi, 4)
    ]

    if any(isequal(oper, o) for o in inv_opers) && isempty(Symbolics.get_variables(x))
        val = eval(Symbolics.toexpr(x))
        for i in eachindex(inv_exacts)
            exact_val = eval(Symbolics.toexpr(inv_exacts[i]))
            if isapprox(exact_val, val, atol=1e-6)
                return inv_exacts[i]
            elseif isapprox(-exact_val, val, atol=1e-6)
                return -inv_exacts[i]
            end
        end
    end

    # add [sin, cos, tan] simplifications in the future?
    return x
end
