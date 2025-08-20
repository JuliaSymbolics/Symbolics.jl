# import Symbolics: symbolic_solve

"""
Returns evolution matrix e^(tD)
"""
evo_mat(D::Matrix{<:Number}, t::Num) = diagm(exp.(t .* diag(D)))

"""
    solve_linear_system(A::Matrix{<:Number}, x0::Vector{<:Number}, t::Num)
Solve linear continuous dynamical system of differential equations of the form Ax = x' with initial condition x0

# Arguments
- `A`: matrix of coefficients
- `x0`: initial conditions vector
- `t`: independent variable

# Returns
vector of symbolic solutions

# Examples
!!! note uses method `symbolic_solve`, so packages `Nemo` and `Groebner` are often required
```jldoctest
julia> @variables t
1-element Vector{Num}:
 t

julia> solve_linear_system([1 0; 0 -1], [1, -1], t) # requires Nemo
2-element Vector{Num}:
   exp(t)
 -exp(-t)

julia> solve_linear_system([-3 4; -2 3], [7, 2], t) # requires Groebner
2-element Vector{Num}:
 (10//1)*exp(-t) - (3//1)*exp(t)
  (5//1)*exp(-t) - (3//1)*exp(t)
```
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
    v = variables(:v, 1:size(A, 1)) # vector of subscripted variables to represent eigenvector
    
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