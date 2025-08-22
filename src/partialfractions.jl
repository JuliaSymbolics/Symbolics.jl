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

Performs partial fraction decomposition for expressions with linear, repeated, or irreducible quadratic factors in the denominator. Can't currently handle irrational roots.

When leading coefficient of the denominator is not 1, it will be factored out and then put back in at the end, often leading to non-integer coefficients in the result. Will return `nothing` if the expression is not a valid polynomial fraction, or if it has irrational roots.

# Examples

```jldoctest
julia> @variables x
1-element Vector{Num}:
 x

julia> partial_frac_decomposition((3x-1) / (x^2 + x - 6), x)
(1//1) / (-2 + x) + (2//1) / (3 + x)

julia> partial_frac_decomposition((4x^3 + 16x + 7)/(x^2 + 4)^2, x) # repeated irreducible quadratic factor
(4x) / (4 + x^2) + 7 / ((4 + x^2)^2)

julia> partial_frac_decomposition((4x^2 - 22x + 7)/((2x+3)*(x-2)^2), x) # non-one leading coefficient
(-3//1) / ((-2 + x)^2) + (2//1) / ((3//2) + x)
```

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

    leading_coeff = coeff_vector(expand(B), x)[end] # of denominator
    
    # already in partial fraction form
    if length(facs) == 1 && only(facs).multiplicity == 1 && degree(A) <= 1
        return expr
    end

    result = []
    c_idx = 0 # index to keep track of which C subscript to use
    if length(facs) == 1
        fac = only(facs)

        if fac.root === nothing # irreducible quadratic factor
            for i = 1:fac.multiplicity
                push!(result, (variable(:𝒞, c_idx+=1)*x + variable(:𝒞, c_idx+=1))/(fac.expr^i)) # (Ax + B)/(x-r)^i
            end
        else
            append!(result, variables(:𝒞, (c_idx+1):(c_idx+=fac.multiplicity)) ./ fac.expr.^(1:fac.multiplicity)) # C1/(x-r) + C2/(x-2)^2 ...
        end
    else
        for fac in facs
            if fac.root === nothing # irreducible quadratic factor
                for i = 1:fac.multiplicity
                    push!(result, (variable(:𝒞, c_idx+=1)*x + variable(:𝒞, c_idx+=1))/(fac.expr^i)) # (Ax + B)/(x-r)^i
                end
                continue
            end

            # cover up method
            other_facs = filter(f -> !isequal(f, fac), facs)
            
            numerator = rationalize(unwrap(substitute(A / prod((f -> f.expr^f.multiplicity).(other_facs)), Dict(x => fac.root)))) # plug in root to expression without its factor in denominator
            push!(result, numerator / fac.expr^fac.multiplicity)

            if fac.multiplicity > 1
                append!(result, variables(:𝒞, (c_idx+1):(c_idx+=fac.multiplicity-1)) ./ fac.expr.^(1:fac.multiplicity-1)) # C1/(x-r) + C2/(x-2)^2 ...
            end
        end
    end

    # no unknowns, so just return
    if isequal(get_variables(sum(result)), [x])
        return sum(result ./ leading_coeff)
    end
        
    lhs = numerator(expr)
    rhs = expand(sum(simplify.(numerator.(result) .* ((B/leading_coeff) ./ denominator.(result))))) # multiply each numerator by the common denominator/its denominator, and sum to get numerator of whole expression

    solution = solve_interms_ofvar(lhs - rhs, x)[1] # solve for unknowns (C's) by looking at coefficients of the polynomial

    # single unknown
    if !(solution isa Dict)
        solution = Dict(variable(:𝒞, 1) => solution)
    end
    
    return sum(substitute.(result, Ref(solution)) ./ leading_coeff) # substitute solutions back in and sum
end

# increasing from 0 to degree n. doesn't skip powers of x like polynomial_coeffs
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
function factorize(expr, x)
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