"""
Represents a linear ordinary differential equation of the form:

dⁿx/dtⁿ + pₙ(t)(dⁿ⁻¹x/dtⁿ⁻¹) + ... + p₂(t)(dx/dt) + p₁(t)x = q(t)

# Fields
- `x`: dependent variable
- `t`: independent variable
- `p`: coefficient functions of `t` ordered in increasing order (p₁, p₂, ...)
- `q`: right hand side function of `t`, without any `x`

# Examples
```jldoctest
julia> using Symbolics

julia> @variables x, t
2-element Vector{Num}:
 x
 t

julia> eq = LinearODE(x, t, [1, 2, 3], 3exp(4t))
(Dt^3)x + (3)(Dt^2)x + (2)(Dt^1)x + (1)(Dt^0)x ~ 3exp(4t)
```
"""
struct LinearODE
    x::Num
    t::Num
    p::AbstractArray
    q::Any
    C::Vector{Num}

    LinearODE(x, t, p, q) = new(x, t, p, q, variables(:C, 1:length(p)))

    function LinearODE(expr, x, t)
        if expr isa Equation
            expr = expr.lhs - expr.rhs
        end

        expr = expand(simplify(expr))

        @assert is_linear_ode(expr, x, t) "Equation must be linear in $x and $t"

        n = _get_der_order(expr, x, t)

        ys = variables(:𝓎, 1:n)
        A, b, islinear = linear_expansion(reduce_order(expr, x, t, ys), ys)

        p = expand.(simplify.(-A[end, 1:end]))
        q = b[end]

        new(x, t, p, q)
    end
end

function is_linear_ode(expr, x, t)
    Dt = Differential(t)
    ys = variables(:𝓎, 1:_get_der_order(expr, x, t))
    n = _get_der_order(expr, x, t)
    @assert n >= 1 "ODE must have at least one derivative"
    
    y_sub = Dict([[(Dt^i)(x) => ys[i+1] for i=0:n-1]; (Dt^n)(x) => variable(:𝒴)])
    expr = substitute(expr, y_sub)

    # isolate (Dt^n)(x)
    f = symbolic_linear_solve(expr, variable(:𝒴), check=false)

    # couldn't isolate
    if f === nothing
        return false
    end

    f = f[1]
    system = [ys[2:n]; f]

    A, b, islinear = linear_expansion(system, ys)
    return islinear && all(isempty.(get_variables.(A, x)))
end

# recursively find highest derivative order in `expr`
function _get_der_order(expr, x, t)
    if !hasderiv(unwrap(expr))
        return 0
    end

    if length(terms(expr)) > 1
        return maximum(_get_der_order.(terms(expr), Ref(x), Ref(t)))
    end

    if length(factors(expr)) > 1
        return maximum(_get_der_order.(factors(expr), Ref(x), Ref(t)))
    end

    return _get_der_order(substitute(expr, Dict(Differential(t)(x) => x)), x, t) + 1
end

Dt(eq::LinearODE) = Differential(eq.t)
order(eq::LinearODE) = length(eq.p)

"""Generates symbolic expression to represent `LinearODE`"""
function get_expression(eq::LinearODE)
    (Dt(eq)^order(eq))(eq.x) + sum([(eq.p[n]) * (Dt(eq)^(n - 1))(eq.x) for n in 1:length(eq.p)]) ~ eq.q
end

function Base.string(eq::LinearODE)
    "(D$(eq.t)^$(order(eq)))$(eq.x) + " *
    join(
        ["($(eq.p[length(eq.p)-n]))(D$(eq.t)^$(length(eq.p)-n-1))$(eq.x)"
         for n in 0:(order(eq) - 1)],
        " + ") * " ~ $(eq.q)"
end

Base.print(io::IO, eq::LinearODE) = print(io, string(eq))
Base.show(io::IO, eq::LinearODE) = print(io, eq)
Base.isequal(eq1::LinearODE, eq2::LinearODE) =
    isequal(eq1.x, eq2.x) && isequal(eq1.t, eq2.t) &&
    isequal(eq1.p, eq2.p) && isequal(eq1.q, eq2.q)

"""Returns true if q(t) = 0 for linear ODE `eq`"""
is_homogeneous(eq::LinearODE) = isempty(Symbolics.get_variables(eq.q))
"""Returns true if all coefficient functions p(t) of `eq` are constant"""
has_const_coeffs(eq::LinearODE) = all(isempty.(Symbolics.get_variables.(eq.p)))
"""Returns homgeneous version of `eq` where q(t) = 0"""
to_homogeneous(eq::LinearODE) = LinearODE(eq.x, eq.t, eq.p, 0)

"""
Returns the characteristic polynomial p of `eq` (must have constant coefficients) in terms of variable `r`

p(D) = Dⁿ + aₙ₋₁Dⁿ⁻¹ + ... + a₁D + a₀I
"""
function characteristic_polynomial(eq::LinearODE, r)
    poly = 0
    @assert has_const_coeffs(eq) "ODE must have constant coefficients to generate characteristic polynomial"
    p = [eq.p; 1] # add implied coefficient of 1 to highest order
    for i in eachindex(p)
        poly += p[i] * r^(i - 1)
    end

    return poly
end

"""
    symbolic_solve_ode(eq::LinearODE)
Symbolically solve a linear ordinary differential equation

# Arguments
- eq: a `LinearODE` to solve

# Returns
Symbolic solution to the ODE

# Supported Methods
- first-order integrating factor
- constant coefficient homogeneous solutions (can handle repeated and complex characteristic roots)
- exponential and resonant response formula particular solutions (for any linear combination of `exp`, `sin`, `cos`, or `exp` times `sin` or `cos` (e.g. `e^2t * cos(-t) + e^-3t + sin(5t))`)
- method of undetermined coefficients particular solutions
- linear combinations of above particular solutions

# Examples

```jldoctest
julia> using Symbolics; import Nemo, SymPy

julia> @variables x, t
2-element Vector{Num}:
 x
 t

# Integrating Factor (note that SymPy is required for integration)
julia> symbolic_solve_ode(LinearODE(x, t, [5/t], 7t))
(C₁ + t^7) / (t^5)

# Constant Coefficients and RRF (note that Nemo is required to find characteristic roots)
julia> symbolic_solve_ode(LinearODE(x, t, [9, -6], 4exp(3t)))
C₁*exp(3t) + C₂*t*exp(3t) + (2//1)*(t^2)*exp(3t)

julia> symbolic_solve_ode(LinearODE(x, t, [6, 5], 2exp(-t)*cos(t)))
C₁*exp(-2t) + C₂*exp(-3t) + (1//5)*cos(t)*exp(-t) + (3//5)*exp(-t)*sin(t)

# Method of Undetermined Coefficients
julia> symbolic_solve_ode(LinearODE(x, t, [-3, 2], 2t - 5))
(11//9) - (2//3)*t + C₁*exp(t) + C₂*exp(-3t)
```
"""
function symbolic_solve_ode(eq::LinearODE)
    homogeneous_solutions = find_homogeneous_solutions(eq)
    
    if is_homogeneous(eq) && homogeneous_solutions !== nothing
        return homogeneous_solutions
    end
    
    particular_solution = find_particular_solution(eq)
    if homogeneous_solutions !== nothing && particular_solution !== nothing
        return homogeneous_solutions + particular_solution
    end
    
    if order(eq) == 1
        return integrating_factor_solve(eq)
    end
end

"""
    symbolic_solve_ode(expr::Equation, x, t)
Symbolically solve an ODE

# Arguments
- expr: a symbolic ODE
- x: dependent variable
- t: independent variable

# Supported Methods
- all methods of solving linear ODEs mentioend for `symbolic_solve_ode(eq::LinearODE)`
- Clairaut's equation
- Bernoulli equations

# Examples

```jldoctest
julia> using Symbolics; import Nemo

julia> @variables x, t
2-element Vector{Num}:
 x
 t

julia> Dt = Differential(t)
Differential(t)

# LinearODE (via constant coefficients and RRF)
julia> symbolic_solve_ode(9t*x - 6*Dt(x) ~ 4exp(3t), x, t)
C₁*exp(3t) + C₂*t*exp(3t) + (2//1)*(t^2)*exp(3t)

# Clairaut's equation
julia> symbolic_solve_ode(x ~ Dt(x)*t - ((Dt(x))^3), x, t)
C₁*t - (C₁^3)

# Bernoulli equations
julia> symbolic_solve_ode(Dt(x) + (4//t)*x ~ t^3 * x^2, x, t)
1 / (C₁*(t^4) - (t^4)*log(t))
```
"""
function symbolic_solve_ode(expr::Equation, x, t)
    clairaut = solve_clairaut(expr, x, t)
    if clairaut !== nothing
        return clairaut
    end

    bernoulli = solve_bernoulli(expr, x, t)
    if bernoulli !== nothing
        return bernoulli
    end

    if is_linear_ode(expr, x, t)
        eq = LinearODE(expr, x, t)
        return symbolic_solve_ode(eq)
    end
end

"""
Find homogeneous solutions of linear ODE `eq` with integration constants of `eq.C`

Currently only works for constant coefficient ODEs
"""
function find_homogeneous_solutions(eq::LinearODE)
    if has_const_coeffs(eq)
        return const_coeff_solve(to_homogeneous(eq))
    end
end

"""
Find a particular solution to linear ODE `eq`

Currently works for any linear combination of exponentials, sin, cos, or an exponential times sin or cos (e.g. e^2t * cos(-t) + e^-3t + sin(5t))
"""
function find_particular_solution(eq::LinearODE)
    # if q has multiple terms, find a particular solution for each and sum together
    terms = Symbolics.terms(eq.q)
    if length(terms) != 1
        solutions = find_particular_solution.(LinearODE.(Ref(eq.x), Ref(eq.t), Ref(eq.p), terms))
        if any(s -> s === nothing, solutions)
            return nothing
        end
        return sum(solutions)
    end

    if has_const_coeffs(eq)
        rrf = resonant_response_formula(eq)
        if rrf !== nothing
            return rrf
        end
        rrf_trig = exp_trig_particular_solution(eq)
        if rrf_trig !== nothing
            return rrf_trig
        end
    end

    undetermined_coeff = method_of_undetermined_coefficients(eq)
    if undetermined_coeff !== nothing
        return undetermined_coeff
    end
end

"""
Returns homogeneous solutions to linear ODE `eq` with constant coefficients

xₕ(t) = C₁e^(r₁t) + C₂e^(r₂t) + ... + Cₙe^(rₙt)
"""
function const_coeff_solve(eq::LinearODE)
    @variables 𝓇
    p = characteristic_polynomial(eq, 𝓇)
    roots = symbolic_solve(p, 𝓇, dropmultiplicity = false)

    # Handle complex + repeated roots
    solutions = exp.(roots * eq.t)
    for i in eachindex(solutions)[1:(end - 1)]
        j = i + 1

        if imag(roots[i]) != 0 && roots[i] == conj(roots[j])
            solutions[i] = exp(real(roots[i] * eq.t)) * cos(imag(roots[i] * eq.t))
            solutions[j] = exp(real(roots[i] * eq.t)) * sin(imag(roots[i] * eq.t))
        end

        while j <= length(solutions) && isequal(roots[i], roots[j])
            solutions[j] *= eq.t # multiply by t for each repetition
            j += 1
        end
    end

    solution = sum(eq.C .* solutions)
    if solution isa Complex && isequal(imag(solution), 0)
        solution = real(solution)
    end

    return solution
end

"""
Solve almost any first order ODE using an integrating factor
"""
function integrating_factor_solve(eq::LinearODE)
    p = eq.p[1] # only p
    v = 0 # integrating factor
    if isempty(Symbolics.get_variables(p))
        v = exp(p * eq.t)
    else
        v = exp(sympy_integrate(p, eq.t))
    end
    solution = (1 / v) * ((isequal(eq.q, 0) ? 0 : sympy_integrate(eq.q * v, eq.t)) + eq.C[1])

    if !isempty(Symbolics.get_variables(solution, variable(:Integral)))
        return nothing
    end
    return expand(Symbolics.sympy_simplify(solution))
end

"""
Returns a, r from q(t)=a*e^(rt) if it is of that form. If not, returns `nothing`
"""
function get_rrf_coeff(q, t)
    facs = _true_factors(q)

    # handle complex r
    # very convoluted, could probably be improved (possibly by making heavier use of @rule)

    # Description of process:
    # only one factor of c*e^((a + bi)t) -> c*cos(bt)e^at + i*c*sin(bt)e^(at)
    # real(factor) / imag(factor) = cos(bt)/sin(bt) - can extract imaginary part b from this
    # then, divide real(factor) = c*cos(bt)e^at by cos(bt) to get c*e^at
    # call self to get c and a, then add back in b
    get_b = Symbolics.Chain([
        (@rule cos(t) / sin(t) => 1), (@rule cos(~b * t) / sin(~b * t) => ~b)])
    if length(facs) == 1 && !isequal(imag(facs[1]), 0) &&
       !isequal(get_b(real(facs[1]) / imag(facs[1])), real(facs[1]) / imag(facs[1]))
        r_im = get_b(real(facs[1]) / imag(facs[1]))
        real_q = real(facs[1]) / cos(r_im * t)
        if isempty(Symbolics.get_variables(real_q, [t]))
            return real_q, r_im * im
        end
        a, r_re = get_rrf_coeff(real(facs[1]) / cos(r_im * t), t)
        return a, r_re + r_im * im
    end

    a = prod(filter(fac -> isempty(Symbolics.get_variables(fac, [t])), facs))

    not_a = filter(fac -> !isempty(Symbolics.get_variables(fac, [t])), facs) # should just be e^(rt)
    if length(not_a) != 1
        return nothing
    end

    der = expand_derivatives(Differential(t)(not_a[1]))
    r = simplify(der / not_a[1])
    if !isempty(Symbolics.get_variables(r, [t]))
        return nothing
    end

    return a, r
end

function _parse_trig(expr, t)
    parse_sin = Symbolics.Chain([(@rule sin(t) => 1), (@rule sin(~x * t) => ~x)])
    parse_cos = Symbolics.Chain([(@rule cos(t) => 1), (@rule cos(~x * t) => ~x)])

    if !isequal(parse_sin(expr), expr)
        return parse_sin(expr), true
    end

    if !isequal(parse_cos(expr), expr)
        return parse_cos(expr), false
    end

    return nothing
end

"""
For finding particular solution when q(t) = a*e^(rt)*cos(bt) (or sin(bt))
"""
function exp_trig_particular_solution(eq::LinearODE)
    facs = _true_factors(eq.q)

    a = prod(filter(fac -> isempty(Symbolics.get_variables(fac, [eq.t])), facs))

    not_a = filter(fac -> !isempty(Symbolics.get_variables(fac, [eq.t])), facs)

    r = nothing
    b = nothing
    is_sin = false

    if length(not_a) == 1 && _parse_trig(not_a[1], eq.t) !== nothing
        r = 0
        b, is_sin = _parse_trig(not_a[1], eq.t)
    elseif length(not_a) != 2
        return nothing
    elseif get_rrf_coeff(not_a[1], eq.t) !== nothing &&
           _parse_trig(not_a[2], eq.t) !== nothing
        r = get_rrf_coeff(not_a[1], eq.t)[2]
        b, is_sin = _parse_trig(not_a[2], eq.t)
    elseif get_rrf_coeff(not_a[2], eq.t) !== nothing &&
           _parse_trig(not_a[1], eq.t) !== nothing
        r = get_rrf_coeff(not_a[2], eq.t)[2]
        b, is_sin = _parse_trig(not_a[1], eq.t)
    else
        return nothing
    end

    # do complex rrf
    # figure out how many times p needs to be differentiated before denominator isn't 0
    k = 0
    @variables 𝓈
    p = characteristic_polynomial(eq, 𝓈)
    Ds = Differential(𝓈)
    while isequal(substitute(expand_derivatives((Ds^k)(p)), Dict(𝓈 => r+b*im)), 0)
        k += 1
    end

    rrf = expand(simplify(a * exp((r + b * im) * eq.t) * eq.t^k /
                           (substitute(expand_derivatives((Ds^k)(p)), Dict(𝓈 => r+b*im)))))

    return is_sin ? imag(rrf) : real(rrf)
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
    @variables 𝓈
    p = characteristic_polynomial(eq, 𝓈)
    Ds = Differential(𝓈)
    while isequal(substitute(expand_derivatives((Ds^k)(p)), Dict(𝓈 => r)), 0)
        k += 1
    end

    return expand(simplify(a * exp(r * eq.t) * eq.t^k /
                           (substitute(expand_derivatives((Ds^k)(p)), Dict(𝓈 => r)))))
end

function method_of_undetermined_coefficients(eq::LinearODE)
    # constant
    p = eq.p[1]
    if isempty(Symbolics.get_variables(p, eq.t)) && isempty(Symbolics.get_variables(eq.q, eq.t))
        return eq.q // p
    end

    # polynomial
    degree = Symbolics.degree(eq.q, eq.t) # just a starting point
    a = Symbolics.variables(:𝒶, 1:degree+1)
    form = sum(a[n]*eq.t^(n-1) for n = 1:degree+1)
    eq_subbed = substitute(get_expression(eq), Dict(eq.x => form))
    eq_subbed = eq_subbed.lhs - eq_subbed.rhs
    eq_subbed = expand_derivatives(eq_subbed)
    
    try
        coeff_solution = solve_interms_ofvar(eq_subbed, eq.t)
    catch
        coeff_solution = nothing
    end
    
    if degree > 0 && coeff_solution !== nothing && !isempty(coeff_solution) && isequal(expand(substitute(eq_subbed, coeff_solution[1])), 0)
        return substitute(form, coeff_solution[1])
    end

    # exponential
    @variables 𝒶
    coeff = get_rrf_coeff(eq.q, eq.t)
    if coeff !== nothing
        r = coeff[2]
        form = 𝒶*exp(r*eq.t)
        eq_subbed = substitute(get_expression(eq), Dict(eq.x => form))
        eq_subbed = expand_derivatives(eq_subbed)
        coeff_solution = symbolic_solve(eq_subbed, 𝒶)
        
        if coeff_solution !== nothing && !isempty(coeff_solution)
            return substitute(form, coeff_solution[1])
        end
    end

    # sin and cos
    # this is a hacky way of doing things
    @variables 𝒶, 𝒷
    @variables 𝒸𝓈, 𝓈𝓃
    parsed = _parse_trig(_true_factors(eq.q)[end], eq.t)
    if parsed !== nothing
        ω = parsed[1]
        form = 𝒶*cos(ω*eq.t) + 𝒷*sin(ω*eq.t)
        eq_subbed = substitute(get_expression(eq), Dict(eq.x => form))
        eq_subbed = expand_derivatives(eq_subbed)
        eq_subbed = expand(substitute(eq_subbed.lhs - eq_subbed.rhs, Dict(cos(ω*eq.t)=>𝒸𝓈, sin(ω*eq.t)=>𝓈𝓃)))
        cos_eq = simplify(sum(filter(term -> !isempty(Symbolics.get_variables(term, 𝒸𝓈)), terms(eq_subbed)))/𝒸𝓈)
        sin_eq = simplify(sum(filter(term -> !isempty(Symbolics.get_variables(term, 𝓈𝓃)), terms(eq_subbed)))/𝓈𝓃)
        if !isempty(Symbolics.get_variables(cos_eq, [eq.t,𝓈𝓃,𝒸𝓈])) || !isempty(Symbolics.get_variables(sin_eq, [eq.t,𝓈𝓃,𝒸𝓈]))
            coeff_solution = nothing
        else
            coeff_solution = symbolic_solve([cos_eq, sin_eq], [𝒶,𝒷])
        end
        
        if coeff_solution !== nothing && !isempty(coeff_solution)
            return substitute(form, coeff_solution[1])
        end
    end
end

function is_solution(solution, eq)
    if solution === nothing
        return false
    end

    expr = substitute(get_expression(eq), Dict(eq.x => solution))
    expr = expand(expand_derivatives(expr.lhs - expr.rhs))
    return isequal(expr, 0)
end

"""
Initial value problem (IVP) for a linear ODE
"""
struct IVP
    eq::LinearODE
    initial_conditions::Vector{Num} # values at t = 0 of nth derivative of x

    function IVP(eq::LinearODE, initial_conditions::Vector{<:Number})
        @assert length(initial_conditions) == order(eq) "# of Initial conditions must match order of ODE"
        new(eq, initial_conditions)
    end
end


function solve_IVP(ivp::IVP)
    general_solution = symbolic_solve_ode(ivp.eq)
    if general_solution === nothing
        return nothing
    end

    eqs = []
    for i in eachindex(ivp.initial_conditions)
        eq::Num = expand_derivatives((Dt(ivp.eq)^(i-1))(general_solution)) - ivp.initial_conditions[i]

        eq = substitute(eq, Dict(ivp.eq.t => 0), fold=false)
        
        # make sure exp, sin, and cos don't evaluate to floats
        exp0 = substitute(exp(ivp.eq.t), Dict(ivp.eq.t => 0), fold=false)
        sin0 = substitute(sin(ivp.eq.t), Dict(ivp.eq.t => 0), fold=false)
        cos0 = substitute(cos(ivp.eq.t), Dict(ivp.eq.t => 0), fold=false)

        eq = expand(simplify(substitute(eq, Dict(exp0 => 1, sin0 => 0, cos0 => 1), fold=false)))
        push!(eqs, eq)
    end

    return expand(simplify(substitute(general_solution, symbolic_solve(eqs, ivp.eq.C)[1])))
end

"""
Solve Clairaut's equation of the form x = x'*t + f(x').

Returns solution of the form x = C*t + f(C) where C is a constant.
"""
function solve_clairaut(expr, x, t)
    Dt = Differential(t)
    rhs = 0
    if isequal(expr.rhs, x)
        rhs = expr.lhs
    elseif isequal(expr.lhs, x)
        rhs = expr.rhs
    else
        return nothing
    end

    terms = Symbolics.terms(rhs)
    matched = false # if expr contains term Dt(x)*t
    f = 0
    for term in terms
        if isequal(term, Dt(x)*t)
            matched = true
        elseif !isempty(Symbolics.get_variables(term, [t]))
            return nothing
        else
            f += term
        end
    end

    if !matched
        return nothing
    end

    C = Symbolics.variable(:C, 1) # constant of integration
    f = substitute(f, Dict(Dt(x) => C))
    if !isempty(Symbolics.get_variables(f, [x]))
        return nothing
    end

    return C*t + f
end

"""
Linearize a Bernoulli equation of the form dx/dt + p(t)x = q(t)x^n into a `LinearODE` of the form dv/dt + (1-n)p(t)v = (1-n)q(t) where v = x^(1-n)
"""
function linearize_bernoulli(expr, x, t, v)
    Dt = Differential(t)

    if expr isa Equation
        expr = expr.lhs - expr.rhs
    end

    terms = Symbolics.terms(expr)

    p = 0
    q = 0
    n = 0
    leading_coeff = 1
    for term in terms
        if Symbolics.hasderiv(Symbolics.value(term))
            facs = _true_factors(term)
            leading_coeff = prod(filter(fac -> !Symbolics.hasderiv(Symbolics.value(fac)), facs))
            if _get_der_order(term//leading_coeff, x, t) != 1
                return nothing
            end
        elseif !isempty(Symbolics.get_variables(term, [x]))
            facs = _true_factors(term)
            x_fac = filter(fac -> !isempty(Symbolics.get_variables(fac, [x])), facs)
            if length(x_fac) != 1
                return nothing
            end

            if isequal(x_fac[1], x)
                p = prod(filter(fac -> isempty(Symbolics.get_variables(fac, [x])), facs))
            else
                n = degree(x_fac[1])
                q = -prod(filter(fac -> isempty(Symbolics.get_variables(fac, [x])), facs))
            end
        end
    end
    
    p //= leading_coeff
    q //= leading_coeff
    
    return LinearODE(v, t, [p*(1-n)], q*(1-n)), n
end

"""
Solve Bernoulli equations of the form dx/dt + p(t)x = q(t)x^n
"""
function solve_bernoulli(expr, x, t)
    @variables 𝓋
    linearized = linearize_bernoulli(expr, x, t, 𝓋)
    if linearized === nothing
        return nothing
    end

    eq, n = linearized

    solution = symbolic_solve_ode(eq)
    if solution === nothing
        return nothing
    end

    return simplify(solution^(1//(1-n)))
end

"""
Solve Bernoulli equations of the form dx/dt + p(t)x = q(t)x^n with initial condition x(0) = x0
"""
function solve_bernoulli(expr, x, t, x0)
    @variables 𝓋
    eq, n = linearize_bernoulli(expr, x, t, 𝓋)

    v0 = x0^(1-n) # convert initial condition from x(0) to v(0)

    ivp = IVP(eq, [v0])
    solution = solve_IVP(ivp)
    if solution === nothing
        return nothing
    end

    return symbolic_solve(solution ~ x^(1-n), x)
end

# takes into account fractions
function _true_factors(expr)
    facs = factors(expr)
    true_facs::Vector{Number} = []
    frac_rule = @rule (~x)/(~y) => [~x, 1/~y]
    for fac in facs
        frac = frac_rule(fac)
        if frac !== nothing && !isequal(frac[1], 1)
            append!(true_facs, _true_factors(frac[1]))
            append!(true_facs, _true_factors(frac[2]))
        else
            push!(true_facs, fac)
        end
    end

    return convert(Vector{Num}, true_facs)
end

"""
    reduce_order(eq, x, t, ys)

Reduce order of an ODE by substituting variables for derivatives to form a system of first order ODEs
"""
function reduce_order(eq, x, t, ys)
    Dt = Differential(t)
    n = _get_der_order(eq, x, t)
    @assert n >= 1 "ODE must have at least one derivative"
    
    # reduction of order
    y_sub = Dict([[(Dt^i)(x) => ys[i+1] for i=0:n-1]; (Dt^n)(x) => variable(:𝒴)])
    eq = substitute(eq, y_sub)
    
    # isolate (Dt^n)(x)
    f = symbolic_linear_solve(eq, variable(:𝒴), check=false)
    @assert f !== nothing "Failed to isolate highest order derivative term"
    f = f[1]
    system = [ys[2:n]; f]

    return system
end

function unreduce_order(expr, x, t, ys)
    Dt = Differential(t)
    rev_y_sub = Dict(ys[i] => (Dt^(i-1))(x) for i in 1:length(ys))

    return substitute(expr, rev_y_sub)
end