# used to represent linear or irreducible quadratic factors
struct Factor
    expr
    root
    multiplicity
    x
end

function Factor(expr, multiplicity, x)
    fac_rule = @rule x + ~r => -~r
    return Factor(expr, fac_rule(expr), multiplicity, x)
end

function Base.isequal(a::Factor, b::Factor)
    return isequal(a.expr, b.expr) && a.multiplicity == b.multiplicity
end

# https://math.mit.edu/~hrm/18.031/pf-coverup.pdf
"""
    partial_frac_decomposition(expr, x)

Performs partial fraction decomposition for expressions with linear, reapeated, or irreducible quadratic factors in the denominator. Can't handle irrational roots or non-one leading coefficients

!!! note that irreducible quadratic and repeated linear factors require the `Groebner` package to solve a system of equations
"""
function partial_frac_decomposition(expr, x)
    A, B = numerator(expr), denominator(expr)

    # check if both numerator and denominator are polynomials
    if !isequal(polynomial_coeffs(A, [x])[2], 0) || !isequal(polynomial_coeffs(B, [x])[2], 0)
        return nothing
    end

    if degree(A) >= degree(B)
        return nothing
    end
    
    facs = factorize(B, x)
    if facs === nothing
        return nothing
    end

    leading_coeff = coeff_vector(expand(B), x)[end] #simplify(B / prod((f -> f.expr^f.multiplicity).(facs)))
    
    if length(facs) == 1 && only(facs).multiplicity == 1 && degree(A) <= 1
        return expr
    end

    result = []
    c_idx = 0
    if length(facs) == 1
        fac = only(facs)
        if fac.root === nothing
            for i = 1:fac.multiplicity
                push!(result, (variable(:C, c_idx+=1)*x + variable(:C, c_idx+=1))/(fac.expr^i))
            end
        else
            append!(result, variables(:C, (c_idx+1):(c_idx+=fac.multiplicity)) ./ fac.expr.^(1:fac.multiplicity))
        end
    else
        for fac in facs
            if fac.root === nothing
                for i = 1:fac.multiplicity
                    push!(result, (variable(:C, c_idx+=1)*x + variable(:C, c_idx+=1))/(fac.expr^i))
                end
                continue
            end

            other_facs = filter(f -> !isequal(f, fac), facs)
            
            numerator = rationalize(unwrap(substitute(A / prod((f -> f.expr^f.multiplicity).(other_facs)), Dict(x => fac.root))))
            push!(result, numerator / fac.expr^fac.multiplicity)

            if fac.multiplicity > 1
                append!(result, variables(:C, (c_idx+1):(c_idx+=fac.multiplicity-1)) ./ fac.expr.^(1:fac.multiplicity-1))
            end
        end
    end
    result
    if isequal(get_variables(sum(result)), [x])
        return sum(result ./ leading_coeff)
    end
        
    lhs::Vector{Rational} = coeff_vector(numerator(expr), x)
    rhs = coeff_vector(expand(sum(simplify.(numerator.(result) .* ((B/leading_coeff) ./ denominator.(result))))), x)
    # rhs = 
    coeff_vector(numerator(simplify(sum(result))), x)

    if length(lhs) > length(rhs)
        rhs = [rhs; zeros(length(lhs)-length(rhs))]
    elseif length(rhs) > length(lhs)
        lhs = [lhs; zeros(length(rhs)-length(lhs))]
    end

    eqs = []
    for i = 1:length(lhs)
        push!(eqs, lhs[i] ~ rhs[i])
    end
    solution = symbolic_solve(eqs, Symbolics.variables(:C, 1:c_idx))[1]

    if !(solution isa Dict)
        solution = Dict(variable(:C, 1) => solution)
    end
    substitute.(result, Ref(solution))
    return sum(substitute.(result, Ref(solution)) ./ leading_coeff)
end

# increasing from 0 to degree n
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

# increasing from 0 to degree of poly
function coeff_vector(poly, x)
    return coeff_vector(poly, x, degree(poly))
end

function count_multiplicities(facs)
    counts = Dict()
    for fac in facs
        if haskey(counts, fac)
            counts[fac] += 1
        else
            counts[fac] = 1
        end
    end

    return counts
end

# for partial fractions, into linear and irreducible quadratic factors
function factorize(expr, x)::Set{Factor}
    roots = symbolic_solve(expr, x, dropmultiplicity=false)

    counts = count_multiplicities(roots)
    facs = Set()

    for root in keys(counts)
        if !isequal(abs(imag(root)), 0)
            fac_expr = expand((x - root)*(x - conj(root)))
            if !isequal(imag(fac_expr), 0)
                @warn "Encountered issue with complex irrational roots"
                return nothing
            end
            push!(facs, Factor(real(fac_expr), counts[root], x))
            continue
        end

        push!(facs, Factor(x - root, root, counts[root], x))
    end

    return facs
end