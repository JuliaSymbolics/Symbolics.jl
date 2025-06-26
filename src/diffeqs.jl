using Groebner, Nemo

"""
Solve first order separable ODE

(mostly me getting used to Symbolics, not super useful in practice)

For example, dx/dt + p(t)x ~ 0
"""
function firstorder_separable_ode_solve(ex, x, t)
    x, t = Symbolics.value(x), Symbolics.value(t)
    p = Symbolics.coeff(ex.lhs, x) # function of t
    P = sympy_integrate(p, t)
    @variables C
    return simplify(C * exp(-P))
end

"""
Returns evolution matrix e^(tD)
"""
evo_mat(D::Matrix{<:Number}, t::Num) = diagm(exp.(t .* diag(D)))

"""
Solve linear continuous dynamical system of differential equations of the form Ax = x' with initial condition x0

# Arguments
- `A`: matrix of coefficients
- `x0`: intial conditions vector
- `t`: independent variable

# Returns
- vector of symbolic solutions
"""
function solve_linear_system(A::Matrix{<:Number}, x0::Vector{<:Number}, t::Num)
    # Check A is square
    if size(A, 1) != size(A, 2)
        throw(ArgumentError("Matrix A must be square."))
    end

    # Check x0 matches size of A
    if size(A, 1) != length(x0)
        throw(ArgumentError("Initial condition vector x0 must match the size of matrix A."))
    end

    if isdiag(A)
        # If A is diagonal, use uncoupled system solver
        return solve_uncoupled_system(A, x0, t)
    end

    S, D = diagonalize(A)

    return simplify.(S * evo_mat(D, t) * rationalize.(S^-1) * x0)
end

"""
Solve a system of uncoupled ODEs of the form:

    dx/dt = A*x
    
where A is a diagonal matrix.
"""
function solve_uncoupled_system(A::Matrix{<:Number}, x0::Vector{<:Number}, t::Num)
    # Check A is diagonal
    if !isdiag(A)
        throw(ArgumentError("Matrix A must be diagonal."))
    end

    return evo_mat(A, t) * x0
end

"""
Diagonalize matrix A, returning matrix S with eigenvectors as columns and diagonal matrix D with eigenvalues
"""
function diagonalize(A::Matrix{<:Number})::Tuple{Matrix{<:Number},Matrix{<:Number}}
    eigs::Eigen = symbolic_eigen(A)
    S = eigs.vectors
    D = diagm(eigs.values)
    return S, D
end


"""
Replacement for `LinearAlgebra.eigen` function that uses symbolic functions to avoid floating-point inaccuracies
"""
function symbolic_eigen(A::Matrix{<:Number})
    @variables λ # eigenvalue
    v = Symbolics.variables(:v, 1:size(A, 1)) # vector of subscripted variables to represent eigenvector
    
    # find eigenvalues first
    p = det(λ*I - A) ~ 0 # polynomial to solve
    values = symbolic_solve(p, λ) # solve polynomial
    
    # then, find eigenvectors
    S::Matrix{Number} = Matrix(I, size(A, 1), 0) # matrix storing vertical eigenvectors
    
    for value in values
        eqs = (value*I - A) * v# .~ zeros(size(A, 1)) # equations to give eigenvectors
        eqs = substitute(eqs, Dict(v[1] => 1)) # set first element to 1 to constrain solution space

        sol = symbolic_solve(eqs[1:end-1], v[2:end]) # solve all but one equation (because of constraining solutions above)

        # parse multivar solutions into Vector (in order)
        if sol[1] isa Dict
            sol = [sol[1][var] for var in v[2:end]]
        end
        vec::Vector{Number} = prepend!(sol, [1]) # add back the 1 (representing v_1) from substitution
        S = [S vec] # add vec to matrix
    end

    return Eigen(values, S)
end

# tests
@variables x, t

Dt = Differential(t)
ex = Dt(x) + 2 * t * x ~ 0
println(firstorder_separable_ode_solve(ex, x, t))
println()

A = [1 0; 0 -1]
x0 = [1, -1]
println(solve_uncoupled_system(A, x0, t))
println()

# commented out below because currently can't handle complex eigenvalues
# A = [-1 -2; 2 -1]
# x0 = [1, -1]
# println(solve_linear_system(A, x0, t))

# println()
# println(symbolic_eigen([-3 4; -2 3])) # should be [1, -1] and [1 2; 1 1] (or equivalent)
# println()
# println(symbolic_eigen([4 -3; 8 -6])) # should be [-2, 0] and [1 2; 3 4] (or equivalent)
# println()
# println(symbolic_eigen([1 -1 0; 1 2 1; -2 1 -1])) # should be [-1, 1, 2] and [-1 -1 -1; -2 0 1; 7 1 1] (or equivalent)
# println()

println(solve_linear_system([-3 4; -2 3], [7, 2], t))
println()
println(solve_linear_system([4 -3; 8 -6], [7, 2], t))
println()
# println(inv(diagonalize([1 -1 0; 1 2 1; -2 1 -1])[1])) --- shows that inv isn't maintaining symbolics
println(solve_linear_system([1 -1 0; 1 2 1; -2 1 -1], [7, 2, 3], t))
println()
