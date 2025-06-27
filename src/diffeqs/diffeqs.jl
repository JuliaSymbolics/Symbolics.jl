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
- ▢ nonhomogeneous with constant coefficients
- ▢ particular solutions (variation of parameters? undetermined coefficients?)
- ▢ [Differential transform method](https://www.researchgate.net/publication/267767445_A_New_Algorithm_for_Solving_Linear_Ordinary_Differential_Equations)
"""
function symbolic_solve_ode(eq::LinearODE)
    if eq.order == 1
        return integrating_factor_solve(eq)
    end

    if is_homogeneous(eq)
        if has_const_coeffs(eq)
            return const_coeff_solve
        end
    end
end

function const_coeff_solve(eq::LinearODE)
    @variables r
    p = characteristic_polynomial(eq, r)
    roots = symbolic_solve(p, r, dropmultiplicity=false)
    return sum(eq.C .* exp.(roots*eq.t))
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
    return Symbolics.sympy_simplify((1/v) * (sympy_integrate(eq.q*v, eq.t) + eq.C[1]))
end