function get_roots_deg1(expression, x)
    subs, filtered_expr = filter_poly(expression, x)
    coeffs, constant = polynomial_coeffs(filtered_expr , [x])

    @assert isequal(sdegree(coeffs, x), 1) "Expected a polynomial of degree 1 in $x, got $expression"

    m = f_numbers(get(coeffs, x, 0))
    c = f_numbers(get(coeffs, x^0, 0))

    root = -c//m
    root = unwrap(ssubs(root, subs))
    return [root]
end

function get_deg2_with_coeffs(a, b, c)
    a, b, c = f_numbers(a), f_numbers(b), f_numbers(c)

    root1 = (-b + term(ssqrt, (b^2 - 4(a*c))))//2a
    root2 = (-b - term(ssqrt, (b^2 - 4(a*c))))//2a

    return [root1, root2]
end

function get_roots_deg2(expression, x)
    # ax^2 + bx + c = 0
    subs, filtered_expr = filter_poly(expression, x)
    coeffs, constant = polynomial_coeffs(filtered_expr, [x])

    @assert isequal(sdegree(coeffs, x), 2) "Expected a polynomial of degree 2 in $x, got $expression"

    results = (f_numbers(unwrap(ssubs(get(coeffs, x^i, 0), subs))) for i in 2:-1:0)
    a, b, c = results

    root1 = (-b + term(ssqrt, (b^2 - 4(a*c))))//2a
    root2 = (-b - term(ssqrt, (b^2 - 4(a*c))))//2a

    return [root1, root2]
end

function get_roots_deg3(expression, x)
    subs, filtered_expr = filter_poly(expression, x)
    coeffs, constant = polynomial_coeffs(filtered_expr, [x])

    @assert isequal(sdegree(coeffs, x), 3) "Expected a polynomial of degree 3 in $x, got $expression"

    results = (f_numbers(unwrap(ssubs(get(coeffs, x^i, 0), subs))) for i in 3:-1:0)
    a, b, c, d = results

    
    Q = (((3*a*c) - b^2))//(9a^2)
    R = ((9*a*b*c - ((27*(a^2)*d)+2b^3)))//(54a^3)
    
    S = term(scbrt, (R + term(ssqrt, (Q^3+R^2))))
    T = term(scbrt, (R - term(ssqrt, (Q^3+R^2))))

    root1 = S + T - (b//(3*a))
    root2 = -((S+T)//2) - (b//(3*a)) + (im*(term(ssqrt, 3))/2)*(S-T)
    root3 = -((S+T)//2) - (b//(3*a)) - (im*(term(ssqrt, 3))/2)*(S-T)

    return [root1, root2, root3]
end



function get_roots_deg4(expression, x)
    subs, filtered_expr = filter_poly(expression, x)
    coeffs, constant = polynomial_coeffs(filtered_expr, [x])

    @assert isequal(sdegree(coeffs, x), 4) "Expected a polynomial of degree 4 in $x, got $expression"

    results = (f_numbers(unwrap(ssubs(get(coeffs, x^i, 0), subs))) for i in 4:-1:0)
    a, b, c, d, e = results

    p = (8(a*c)-3(b^2))//(8(a^2))

    q = (b^3 - 4(a*b*c) + 8(d*a^2))//(8*a^3)

    r = (-3(b^4) + 256(e*a^3) - 64(d*b*a^2) + 16(c*(b^2)*a))//(256*a^4)

    m = gensym()
    m = (@variables $m)[1]
    eq_m = 8m^3 + 8(p)*m^2 + (2(p^2) - 8r)m - q^2

    roots_m = solve_univar(eq_m, m)
    m = 0

    # Yassin: this thing is a problem for parametric
    for root in roots_m
        vars = get_variables(root)
        if isequal(vars, []) && !isequal(eval(toexpr(root)), 0)
            m = Symbolics.unwrap(copy(Symbolics.wrap(root)))
            break
        end
    end
    if isequal(m, 0)
        @info "Assuming $(roots_m[1] != 0)"
        m = roots_m[1]
    end

    arr = get_yroots(m, p, q)
    for i in eachindex(arr)
        arr[i] -= unwrap(b//4a)
    end

    return arr
end

function get_yroots(m, p, q)
    a = 1
    b1 = term(ssqrt,  2m)
    c1 = (p//2) + m - (q//(2*term(ssqrt, 2m)))
    b2 = -term(ssqrt, 2m)
    c2 = (p//2) + m + (q//(2*term(ssqrt, 2m)))

    root1, root2 = get_deg2_with_coeffs(a, b1, c1)
    root3, root4 = get_deg2_with_coeffs(a, b2, c2)
    return [root1, root2, root3, root4]
end

function get_roots(expression, x)
    @assert is_singleton(unwrap(x)) "Expected a variable, got $x"

    subs, filtered_expr = filter_poly(expression, x)
    coeffs, constant = polynomial_coeffs(filtered_expr, [x])
    @assert isequal(constant, 0) "Expected a polynomial in $x, got $expression"

    degree = sdegree(coeffs, x)

    if degree == 0 && isequal(expression, 0)
        return [x]
    elseif degree == 0 && !isequal(expression, 0)
        return [] # no roots!
    end



    if degree == 1
        return get_roots_deg1(expression, x)
    end

    if degree == 2
        return get_roots_deg2(expression, x)
    end

    if degree == 3
        return get_roots_deg3(expression, x)
    end

    if degree == 4
        return get_roots_deg4(expression, x)
    end

end

