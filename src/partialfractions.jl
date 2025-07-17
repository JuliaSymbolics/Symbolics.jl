# https://dl.acm.org/doi/pdf/10.1145/800204.806314
function partial_frac_decompose(expr)
    frac_rule = @rule ~a/~b => (~a, ~b)
    A, B = frac_rule(expr)
    
    X = []
    Z = []
    b, Y = factor_use_nemo(B)
    Y = expand_factor_list(Y)
    k = length(Y)
    
    if k == 1
        return [A], [B], [1]
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
    
    E = 1 # TODO
    F = 1 # TODO

    G = F*E^(-1) # TODO solve EG = F for G

    w = b * G[0]

    for i = 1:k
        m = degree(B[i])
        
        if m == 0
            push!(X, 0)
            push!(Z, 1)
            continue
        end

        n = i*m
        
    end
end

function expand_factor_list(facs)
    result = []
    expand_rule = @rule ~a^~n => fill(~a, ~n)
    for fac in facs
        expanded = expand_rule(fac)
        if expanded === nothing
            push!(result, fac)
        else
            append!(result, expanded)
        end
    end

    return result
end