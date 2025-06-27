using Symbolics
import Symbolics: value, coeff, sympy_integrate

struct LinearODE
    # dⁿx/dtⁿ + pₙ(t)(dⁿ⁻¹x/dtⁿ⁻¹) + ... + p₂(t)(dx/dt) + p₁(t)x = q(t)

    x::SymbolicUtils.Symbolic # dependent variable
    t::SymbolicUtils.Symbolic # independent variable
    p::AbstractArray # coefficient functions of t ordered in increasing order (p₁, p₂, ...)
    q # right hand side function of t, without any x

    LinearODE(x::Num, t::Num, p, q) = new(value(x), value(t), p, q)
end


is_homogeneous(eq::LinearODE) = isempty(Symbolics.get_variables(eq.q))

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