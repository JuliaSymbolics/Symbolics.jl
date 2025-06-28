using Symbolics
import Symbolics: value, coeff, sympy_integrate

struct LinearODE
    # dⁿx/dtⁿ + pₙ(t)(dⁿ⁻¹x/dtⁿ⁻¹) + ... + p₂(t)(dx/dt) + p₁(t)x = q(t)

    x::SymbolicUtils.Symbolic # dependent variable
    t::SymbolicUtils.Symbolic # independent variable
    p::AbstractArray # coefficient functions of t ordered in increasing order (p₁, p₂, ...)
    q # right hand side function of t, without any x
    Dt::Differential
    order::Int
    C::Vector{Num} # constants

    LinearODE(x::Num, t::Num, p, q) = new(value(x), value(t), p, q, Differential(t), length(p), variables(:C, 1:length(p)))
    LinearODE(x, t, p, q) = new(x, t, p, q, Differential(t), length(p), variables(:C, 1:length(p)))
end

get_expression(eq::LinearODE) = (eq.Dt^eq.order)(eq.x) + sum([(eq.p[n])*(eq.Dt^(n-1))(eq.x) for n = 1:length(eq.p)]) ~ eq.q

Base.print(io::IO, eq::LinearODE) = print(io, "(D$(eq.t)^$(eq.order))$(eq.x) + " * join(["($(eq.p[length(eq.p)-n]))(D$(eq.t)^$(length(eq.p)-n-1))$(eq.x)" for n = 0:(eq.order-1)], " + ") * " ~ $(eq.q)")
Base.show(io::IO, eq::LinearODE) = print(io, eq)

is_homogeneous(eq::LinearODE) = isempty(Symbolics.get_variables(eq.q))
has_const_coeffs(eq::LinearODE) = all(isempty.(Symbolics.get_variables.(eq.p)))

to_homogeneous(eq::LinearODE) = LinearODE(eq.x, eq.t, eq.p, 0)

function characteristic_polynomial(eq::LinearODE, r)
    poly = 0
    @assert has_const_coeffs(eq) "ODE must have constant coefficients to generate characteristic polynomial"
    p = [eq.p; 1] # add implied coefficient of 1 to highest order
    for i in eachindex(p)
        poly += p[i] * r^(i-1)
    end

    return poly
end

"""
Symbolically solve a linear ODE

Cases handled:
- ☑ first order
- ☑ homogeneous with constant coefficients
- ◩ nonhomogeneous with constant coefficients
    - ☑ ERF + RRF
- ▢ particular solutions (variation of parameters? undetermined coefficients?)
- ▢ [Differential transform method](https://www.researchgate.net/publication/267767445_A_New_Algorithm_for_Solving_Linear_Ordinary_Differential_Equations)
"""
function symbolic_solve_ode(eq::LinearODE)
    if eq.order == 1
        return integrating_factor_solve(eq)
    end

    if is_homogeneous(eq)
        if has_const_coeffs(eq)
            return const_coeff_solve(eq)
        end
    end

    if has_const_coeffs(eq)
        rrf = resonant_response_formula(eq)
        if rrf !== nothing
            return const_coeff_solve(to_homogeneous(eq)) + rrf
        end
    end
end

function const_coeff_solve(eq::LinearODE)
    @variables r
    p = characteristic_polynomial(eq, r)
    roots = symbolic_solve(p, r, dropmultiplicity=false)

    # Handle repeated roots
    solutions = exp.(roots*eq.t)
    for i in eachindex(solutions)[1:end-1]
        j = i+1
        while j <= length(solutions) && isequal(solutions[i], solutions[j])
            solutions[j] *= eq.t^(j-i) # multiply by t for each repetition
            j+=1
        end
    end

    return sum(eq.C .* solutions)
end

"""
Solve almost any first order ODE using an integrating factor
"""
function integrating_factor_solve(eq::LinearODE)
    p = eq.p[1] # only p
    v = 0 # integrating factor
    if isempty(Symbolics.get_variables(p))
        v = exp(p*eq.t)
    else
        v = exp(sympy_integrate(p, eq.t))
    end
    return Symbolics.sympy_simplify((1/v) * ((isequal(eq.q, 0) ? 0 : sympy_integrate(eq.q*v, eq.t)) + eq.C[1]))
end

"""
Returns a, r from q(t)=a*e^(rt) if it is of that form. If not, returns `nothing`
"""
function get_rrf_coeff(q, t)
    facs = factors(q)
    
    # handle complex r
    # very convoluted, could probably be improved (possibly by making heavier use of @rule)

    # Description of process:
    # only one factor of c*e^((a + bi)t) -> c*cos(bt)e^at + i*c*sin(bt)e^(at)
    # real(factor) / imag(factor) = cos(bt)/sin(bt) - can extract imaginary part b from this
    # then, divide real(factor) = c*cos(bt)e^at by cos(bt) to get c*e^at
    # call self to get c and a, then add back in b
    get_b = @rule cos(~b*t) / sin(~b*t) => ~b
    if length(facs) == 1 && !isequal(imag(facs[1]), 0) && get_b(real(facs[1])/imag(facs[1])) !== nothing
        r_im = get_b(real(facs[1])/imag(facs[1]))
        real_q = real(facs[1]) / cos(r_im*t)
        if isempty(Symbolics.get_variables(real_q, t))
            return real_q, r_im*im
        end
        a, r_re = get_rrf_coeff(real(facs[1]) / cos(r_im*t), t)
        return a, r_re + r_im*im
    end

    a = prod(filter(fac -> isempty(Symbolics.get_variables(fac, [t])), facs))
    
    not_a = filter(fac -> !isempty(Symbolics.get_variables(fac, [t])), facs) # should just be e^(rt)
    if length(not_a) != 1
        return nothing
    end

    der = expand_derivatives(Differential(t)(not_a[1]))
    r = simplify(der / not_a[1])
    if !isempty(Symbolics.get_variables(r, t))
        return nothing
    end

    return a, r
end

"""
Returns a particular solution to a constant coefficient ODE with q(t) = a*e^(rt)

Exponential Response Formula: x_p(t) = a*e^(rt)/p(r) where p(r) is characteristic polynomial

Resonant Response Formula: If r is a characteristic root, multiply by t and take the derivative of p (possibly multiple times)
"""
function resonant_response_formula(eq::LinearODE)
    @assert has_const_coeffs(eq)

    # get a and r from q = a*e^(rt)
    rrf_coeff = get_rrf_coeff(eq.q, eq.t)
    if rrf_coeff === nothing
        return nothing
    end
    a, r = rrf_coeff
    
    # figure out how many times p needs to be differentiated before denominator isn't 0
    k = 0
    @variables s
    p = characteristic_polynomial(eq, s)
    Ds = Differential(s)
    while isequal(substitute(expand_derivatives((Ds^k)(p)), Dict(s => r)), 0)
        k += 1
    end

    return (eq.q*eq.t^k) / (substitute(expand_derivatives((Ds^k)(p)), Dict(s => r)))
end