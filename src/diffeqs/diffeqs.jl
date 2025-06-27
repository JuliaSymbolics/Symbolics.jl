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

    LinearODE(x::Num, t::Num, p, q) = new(value(x), value(t), p, q, Differential(t), length(p)+1)
    LinearODE(x, t, p, q) = new(x, t, p, q, Differential(t), length(p)+1)
end

get_expression(eq::LinearODE) = (eq.Dt^eq.order)(eq.x) + sum([(eq.p[n])*(eq.Dt^n)(eq.x) for n = 1:length(eq.p)]) ~ eq.q

Base.print(io::IO, eq::LinearODE) = print(io, "(D$(eq.t)^$(eq.order))$(eq.x) + " * join(["($(eq.p[length(eq.p)-n]))(D$(eq.t)^$(length(eq.p)-n))$(eq.x)" for n = 0:(eq.order-2)], " + ") * " ~ $(eq.q)")
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

function homogeneous_solve(eq::LinearODE)
    @variables r
    p = characteristic_polynomial(eq, r)
    roots = symbolic_solve(p, r, dropmultiplicity=false)
    @variables C[1:degree(p, r)]
    return sum(Symbolics.scalarize(C) .* exp.(roots*eq.t))
end

"""
Solve first order separable ODE

(mostly me getting used to Symbolics, not super useful in practice)

For example, dx/dt + p(t)x ~ 0
"""
function firstorder_separable_ode_solve(ex, x, t)
    x, t = Symbolics.value(x), Symbolics.value(t)
    p = Symbolics.coeff(ex.lhs, x) # function of t
    P = Symbolics.sympy_integrate(p, t)
    @variables C
    return simplify(C * exp(-P))
end