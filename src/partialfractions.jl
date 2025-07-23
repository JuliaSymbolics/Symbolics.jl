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
    
    facs = factorize(B, x)
    v = simplify(B / prod((f -> f.expr^f.multiplicity).(facs)))
    
    if length(facs) == 1 && only(facs).multiplicity == 1 && degree(A) <= 1
        return expr
    end

    result = 0

    c_idx = 0
    if length(facs) == 1
        fac = only(facs)
        if fac.root === nothing
            for i = 1:fac.multiplicity
                result += (variable(:C, c_idx+=1)*x + variable(:C, c_idx+=1))/(fac.expr^i)
            end
        else
            result += sum(variables(:C, (c_idx+1):(c_idx+=fac.multiplicity)) ./ fac.expr.^(1:fac.multiplicity))
        end
    else
        for fac in facs
            if fac.root === nothing
                for i = 1:fac.multiplicity
                    result += (variable(:C, c_idx+=1)*x + variable(:C, c_idx+=1))/(fac.expr^i)
                end
                continue
            end

            other_facs = filter(f -> !isequal(f, fac), facs)
            
            numerator = rationalize(unwrap(substitute(A / prod((f -> f.expr^f.multiplicity).(other_facs)), Dict(x => fac.root))))
            result += numerator / fac.expr^fac.multiplicity

            if fac.multiplicity > 1
                result += sum(variables(:C, (c_idx+1):(c_idx+=fac.multiplicity-1)) ./ fac.expr.^(1:fac.multiplicity-1))
            end
        end
    end

    if isequal(get_variables(result), [x])
        return expand(result/v)
    end
        
    lhs::Vector{Rational} = coeff_vector(numerator(expr), x)
    rhs = coeff_vector(numerator(simplify(result)), x)

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
    
    return expand(substitute(result, solution)/v)
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
            fac_expr = real(expand(real((x - root)*(x - conj(root)))))
            push!(facs, Factor(fac_expr, counts[root], x))
            continue
        end

        push!(facs, Factor(x - root, root, counts[root], x))
    end

    return facs
end