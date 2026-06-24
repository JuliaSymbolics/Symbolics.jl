using SymbolicUtils: unsorted_arguments, quick_cancel, Div
# Fix (x+5)/((x-1)^(-2) * (x-2))
function partial_fraction(expr, var)
    expr = unwrap(expr)
    var = unwrap(var)

    expr isa Div || return expr
    den, num = expr.den, expr.num

    den_total_deg = 0
    den_args = unsorted_arguments(den)
    sol = []

    if operation(den) == (*)
        for factor in den_args
            if factor isa Pow
                b, a, islinear = linear_expansion(factor.base, var)
                islinear || return expr
                isinteger(factor.exp) || return expr
                den_total_deg += factor.exp
                push!(sol, (b, -a // b))
            elseif factor isa Add
                b, a, islinear = linear_expansion(factor, var)
                if islinear
                    push!(sol, -a // b)
                    den_total_deg += 1
                else
                    return expr
                end
            end
        end
        den_total_deg <= 1 && return expr
        degree(num) < den_total_deg || return expr

        out = 0
        for (s, factor) in zip(sol, den_args)
            # f = (ax + b) * g
            g = simplify_fractions(factor * expr)
            if factor isa Add
                out += substitute(g, var=>s) / factor
            else # factor isa Pow
                b, s = s
                count = 0
                factorial = 1
                current_exp = factor.exp
                while (current_exp > 0)
                    numer = substitute(g, var=>s)
                    out += numer / (factor.base^current_exp * factorial * b^count)
                    current_exp -= 1
                    count += 1
                    factorial *= count
                    g = simplify_fractions(derivative(g, var))
                end
            end
        end
        return out
    else
        return expr
    end
end
