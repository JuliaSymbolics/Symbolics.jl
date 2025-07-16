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

    args = arguments(x)

    # (X)^0 => 1
    if oper === (^) && isequal(args[2], 0) && !isequal(args[1], 0)
        return 1
    end

    # (X)^1 => X
    if oper === (^) && isequal(args[2], 1)
        return args[1]
    end

    # (0)^X => 0
    if oper === (^) && isequal(args[1], 0) && !isequal(args[2], 0)
        return 0
    end

    # y / 0 => Inf
    if oper === (/) && !isequal(numerator(x), 0) && isequal(denominator(x), 0)
        return Inf
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

    x = convert_consts(x)

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


inv_exacts = [0, Symbolics.term(*, pi),
        Symbolics.term(/,pi,3),
        Symbolics.term(/, pi, 2),
        Symbolics.term(/, Symbolics.term(*, 2, pi), 3),
        Symbolics.term(/, pi, 6),
        Symbolics.term(/, Symbolics.term(*, 5, pi), 6),
        Symbolics.term(/, pi, 4)
]
inv_evald = Symbolics.symbolic_to_float.(inv_exacts)

const inv_pairs = collect(zip(inv_exacts, inv_evald))
"""
    function convert_consts(x)
This function takes BasicSymbolic terms as input (x) and attempts
to simplify these basic symbolic terms using known values.
Currently, this function only supports inverse trigonometric functions.

## Examples
```jldoctest
julia> Symbolics.convert_consts(Symbolics.term(acos, 0))
π / 2

julia> Symbolics.convert_consts(Symbolics.term(atan, 0))
0

julia> Symbolics.convert_consts(Symbolics.term(atan, 1))
π / 4
```
"""
function convert_consts(x)
    !iscall(x) && return x

    oper = operation(x)
    inv_opers = [asin, acos, atan]

    if any(isequal(oper, o) for o in inv_opers) && isempty(Symbolics.get_variables(x))
        val = Symbolics.symbolic_to_float(x)
        for (exact, evald) in inv_pairs
            if isapprox(evald, val)
                return exact
            elseif isapprox(-evald, val)
                return -exact
            end
        end
    end

    # add [sin, cos, tan] simplifications in the future?
    return x
end
