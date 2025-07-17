# https://dl.acm.org/doi/pdf/10.1145/800204.806314
function partial_frac_decomposition(expr, x)
    frac_rule = @rule ~a/~b => (~a, ~b)
    A, B = frac_rule(expr)
    
    X = []
    Z = []
    b = gcd(coeff_vector(B, x)...)
    Y = x .- symbolic_solve(B, x, dropmultiplicity=false)
    k = length(Y)
    
    if k == 1
        return expr
    end

    # idk what this is for
    # for i = 1:k-1
    #     if degree(B[i]) == 0
    #         push!(X, A)
    #         push!(Z, b)
    #         for j = 1:k-1
    #             push!(X, 0)
    #             push!(Z, 1)
    #         end
    #         return X, Y, Z
    #     end
    # end
    
    E = zeros(degree(B), degree(B))
    F = coeff_vector(A, x)

    for i = 1:k
        E[:, i] = coeff_vector(simplify(B/Y[i]), x, degree(B) - 1)
    end
    @show E
    Gbar = rationalize.(inv(E), tol=1e-6)*F
    G0 = Gbar[1]
    G = Gbar[2:end] / G0

    w = b * G0

    j = 1
    for i = 1:k
        m = degree(Y[i])
        
        if m == 0
            push!(X, 0)
            push!(Z, 1)
            continue
        end

        Ai = G[j:j+m] .* x.^[0:m]
        h = gcd(w, coeff_vector(Ai, x)...)
        vi = w / h
        Ai /= h

        if vi < 0
            vi *= -1
            Ai *= -1
        end

        push!(X, Ai)
        push!(Z, vi)

        j += m
    end

    return reverse(X) ./ (Y .* reverse(Z))
end

function coeff_vector(poly, x, n)
    coeff_dict = polynomial_coeffs(poly, [x])[1]
    vec = []
    for i = 0:n
        if x^i in keys(coeff_dict)
            push!(vec, coeff_dict[x^i])
        else
            push!(vec, 0)
        end
    end

    return vec
end

function coeff_vector(poly, x)
    return coeff_vector(poly, x, degree(poly))
end