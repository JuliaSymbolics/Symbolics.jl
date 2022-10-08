function nterms(t)
    if istree(t)
        return reduce(+, map(nterms, arguments(t)), init=0)
    else
        return 1
    end
end

# Soft pivoted
# Note: we call this function with a matrix of Union{SymbolicUtils.Symbolic, Any}
function _sym_lu(A)
    SINGULAR = typemax(Int)
    m, n = size(A)
    F = map(x->x isa RCNum ? x : Num(x), A)
    minmn = min(m, n)
    p = Vector{LinearAlgebra.BlasInt}(undef, minmn)
    lead = 1
    leads = zeros(Int, minmn)
    rank = 0
    for k = 1:minmn
        kp = k # pivot index

        # search for the expression different from zero with the least terms
        amin = SINGULAR
        for j = lead:n
            for i = k:m # search first by columns
                absi = _iszero(F[i, j]) ? SINGULAR : nterms(F[i, j])
                if absi < amin
                    kp = i
                    amin = absi
                end
            end
            # break when pivot found
            if amin != SINGULAR
                lead = j
                leads[k] = lead
                break
            end
        end

        p[k] = kp

        # break from function as the reduced echelon form has been
        # reached, but fill `p`
        if amin == SINGULAR && !(amin isa Symbolic) && (amin isa Number)
            for i = k+1:minmn
                p[i] = i
            end
            break
        end
        rank = k

        # swap rows
        if k != kp
            for j = 1:n
                F[k, j], F[kp, j] = F[kp, j], F[k, j]
            end
        end

        # set values for L matrix
        c = F[k, lead]
        for i = k+1:m
            F[i, k] = F[i, lead] / c
            if lead != k
                F[i, lead] = zero(Num)
            end
        end

        # substract the row from every other, traverse first by colums
        # we start from lead+1 to avoid chaing the leading value on the column
        for j = lead+1:n
            for i = k+1:m
                # this line occupies most of the time, distributed in the
                # following methods
                #   - `*(::Num, ::Num)` dynamic dispatch
                #   - `-(::Num, ::Num)` dynamic dispatch
                F[i, j] = F[i, j] - F[i, k] * F[k, j]
            end
        end

        # advance the lead by one
        lead = lead + 1
    end
    return F, p, filter!(!iszero, leads), rank
end
function sym_lu(A; check=true)
    F, p, leads, rank = _sym_lu(A)
    info = rank == minimum(size(A)) ? 0 : rank
    check && LinearAlgebra.checknonsingular(info)
    LU(F, p, convert(LinearAlgebra.BlasInt, info))
end

# soft pivoted reduced row echelon form for extended matrix
function sym_rref(A, b)
    SINGULAR = typemax(Int)
    m, n = size(A)
    F = map(x->x isa Num ? x : Num(x), A)
    b = map(x->x isa Num ? x : Num(x), b)
    minmn = min(m, n)
    lead = 1
    leads = zeros(Int, minmn)
    rank = 0
    for k = 1:m
        kp = k # pivot index

        # search for the expression different from zero with the least terms
        amin = SINGULAR
        for j = lead:n
            for i = k:m # search first by columns
                absi = _iszero(F[i, j]) ? SINGULAR : nterms(F[i, j])
                if absi < amin
                    kp = i
                    amin = absi
                end
            end
            # break when pivot found
            if amin != SINGULAR
                lead = j
                leads[k] = lead
                break
            end
        end

        # normally, here we would save the permutation in an array
        # but this is not needed as we have the extended matrix

        # break from function as the reduced echelon form has been reached
        if amin == SINGULAR && !(amin isa Symbolic) && (amin isa Number)
            break
        end
        rank = k

        # swap rows, only needed to swap lead:end
        if k != kp
            for j = lead:n
                F[k, j], F[kp, j] = F[kp, j], F[k, j]
            end
        end

        # normalize the current row
        c = F[k, lead]
        for j = lead:n
            F[k, j] = F[k, j] / c
        end
        b[k] = b[k] / c

        # substract the row form every other, traverse first by colums
        for j = lead+1:n
            for i = 1:m
                if i != k
                    # this line occupies most of the time, distributed in the
                    # following methods
                    #   - `*(::Num, ::Num)` dynamic dispatch
                    #   - `-(::Num, ::Num)` dynamic dispatch
                    F[i, j] = F[i, j] - F[i, lead] * F[k, j]
                end
            end
        end
        # substract the row in the extended part of the matrix
        for i = 1:m
            if i != k
                b[i] = b[i] - F[i, lead] * b[k]
            end
        end
        # zero the lead column
        for i = 1:m
            if i != k
                F[i, lead] = zero(Num)
            end
        end

        # advance the lead by one
        lead = lead + 1
    end
    return F, b, filter!(!iszero, leads), rank
end

# convert an upper extended matrix to rref using leads
# modifies `U` and `b`
function _sym_urref!(U, b, leads)
    m, n = size(U)
    for k = 1:length(leads)
        lead = leads[k]

        c = U[k, lead]
        for j = lead:n
            U[k, j] = U[k, j] / c
        end
        b[k] = b[k] / c

        for j = lead+1:n
            for i = 1:k-1
                # this line occupies most of the time, distributed in the
                # following causes
                #   - `*(::Num, ::Num)` dynamic dispatch
                #   - `-(::Num, ::Num)` dynamic dispatch
                U[i, j] = U[i, j] - U[k, j] * U[i, lead]
            end
        end
        for i = 1:m
            if i != k
                b[i] = b[i] - b[k] * U[i, lead]
            end
        end
        for i = 1:m
            if i != k
                U[i, lead] = zero(eltype(U))
            end
        end
    end
    U, b
end

function _factorize(A, b)
    m, n = size(A)
    if m <= n
        F, p, leads, rank = _sym_lu(A)
        return F, b, p, leads, rank
    else
        F, b, leads, rank = sym_rref(A, b)
        return F, b, Int[], leads, rank
    end
end

# Given a vector of equations and a
# list of independent variables,
# return the coefficient matrix `A` and a
# vector of constants (possibly symbolic) `b` such that
# A \ b will solve the equations for the vars
function A_b(eqs::AbstractArray, vars::AbstractArray, check)
    exprs = rhss(eqs) .- lhss(eqs)
    if check
        for ex in exprs
            @assert isaffine(ex, vars)
        end
    end
    A = jacobian(exprs, vars)
    b = A * vars - exprs
    A, b
end

"""
$(TYPEDSIGNATURES)

Solve equation(s) `eqs` for a set of variables `vars`.

Assumes `length(eqs) == length(vars)`

Currently only works if all equations are linear. `check` if the expr is linear
w.r.t `vars`.

# Examples
```julia
julia> @variables x y
2-element Vector{Num}:
 x
 y

julia> Symbolics.solve_for(x + y ~ 0, x)
-y

julia> Symbolics.solve_for([x + y ~ 0, x - y ~ 2], [x, y])
2-element Vector{Float64}:
  1.0
 -1.0
```
"""
function solve_for(eq, var; simplify=false, check=true) # scalar case
    # simplify defaults for `false` as canonicalization should handle most of
    # the cases.
    a, b, islinear = linear_expansion(eq, var)
    check && @assert islinear
    islinear || return nothing
    # a * x + b = 0
    if eq isa AbstractArray && var isa AbstractArray
        x = _solve(a, -b, var, simplify)
    else
        x = a \ -b
    end
    simplify || return x
    if x isa AbstractArray
        SymbolicUtils.simplify.(simplify_fractions.(x))
    else
        SymbolicUtils.simplify(simplify_fractions(x))
    end
end
solve_for(eq::Equation, var::T; x...) where {T<:AbstractArray} = solve_for([eq],var, x...)
solve_for(eq::T, var::Num; x...) where {T<:AbstractArray} = first(solve_for(eq,[var], x...))

const ℝ = (identity)((Symbolics.wrap)((SymbolicUtils.setmetadata)((SymbolicUtils.Sym){Real}(:ℝ), Symbolics.VariableSource, (:variables, :ℝ))))

function _solve(A::AbstractMatrix, b::AbstractArray, vars, do_simplify)
    A = Num.(SymbolicUtils.quick_cancel.(A))
    b = Num.(SymbolicUtils.quick_cancel.(b))
    m, n = size(A)
    minmn = min(m, n)
    F, b, ipiv, leads, rank = _factorize(A, b)
    if m == n && rank == minmn
        info = 0
        sol = value.(LU(F, ipiv, convert(LinearAlgebra.BlasInt, info)) \ b)
    else
        # check for consistency
        for i ∈ rank+1:m
            if !iszero(b[i])
                throw(ArgumentError("Inconsistent linear system"))
            end
        end
        _sol = Vector{Num}(undef, n)
        if m <= n
            # system is of the form Ax = (LU)x = L(Ux) = Lx' = b[p]
            # with L being square `UnitLowerTriangular`

            # first solve Lx' = b[p], -----------------------------------+
            p = LinearAlgebra.ipiv2perm(ipiv, m)                       # |
            L = UnitLowerTriangular(F[1:m, 1:m])                       # |
            _x = symsub!(L, b[p])                                      # |
                                                                       # |
            # then solve Ux = x', first converting [U|x'] to rref, <-----+
            U = triu!(F)                                               # |
            F, b = _sym_urref!(U, _x, leads)                           # |
        else                                                           # |
        end                                                            # |
        # and later filling the values for all the variables <-----------+
        freeidx = setdiff(1:n, leads) # indices for free variables
        for i in freeidx
            _sol[i] = ℝ
        end
        F[CartesianIndex.(1:length(leads), leads)] .= 0
        for (i, v) ∈ enumerate(leads)
            _sol[v] = b[i] - view(F, i, :) ⋅ vars
        end

        sol = value.(_sol)
    end
    do_simplify ? SymbolicUtils.simplify_fractions.(sol) : sol
end

LinearAlgebra.ldiv!(A::UpperTriangular{<:Union{Symbolic,RCNum}}, b::AbstractVector{<:Union{Symbolic,RCNum}}, x::AbstractVector{<:Union{Symbolic,RCNum}} = b) = symsub!(A, b, x)
function symsub!(A::UpperTriangular, b::AbstractVector, x::AbstractVector = b)
    LinearAlgebra.require_one_based_indexing(A, b, x)
    n = size(A, 2)
    if !(n == length(b) == length(x))
        throw(DimensionMismatch("second dimension of left hand side A, $n, length of output x, $(length(x)), and length of right hand side b, $(length(b)), must be equal"))
    end
    @inbounds for j = n:-1:1
        _iszero(A.data[j,j]) && throw(SingularException(j))
        xj = x[j] = b[j] / A.data[j,j]
        for i = j-1:-1:1
            sub = _isone(xj) ? A.data[i,j] : A.data[i,j] * xj
            if !_iszero(sub)
                b[i] -= sub
            end
        end
    end
    x
end

LinearAlgebra.ldiv!(A::UnitLowerTriangular{<:Union{Symbolic,RCNum}}, b::AbstractVector{<:Union{Symbolic,RCNum}}, x::AbstractVector{<:Union{Symbolic,RCNum}} = b) = symsub!(A, b, x)
function symsub!(A::UnitLowerTriangular, b::AbstractVector, x::AbstractVector = b)
    LinearAlgebra.require_one_based_indexing(A, b, x)
    n = size(A, 2)
    if !(n == length(b) == length(x))
        throw(DimensionMismatch("second dimension of left hand side A, $n, length of output x, $(length(x)), and length of right hand side b, $(length(b)), must be equal"))
    end
    @inbounds for j = 1:n
        xj = x[j] = b[j]
        for i = j+1:n
            sub = _isone(xj) ? A.data[i,j] : A.data[i,j] * xj
            if !_iszero(sub)
                b[i] -= sub
            end
        end
    end
    x
end

minor(B, j) = B[2:end, 1:size(B,2) .!= j]
function LinearAlgebra.det(A::AbstractMatrix{<:RCNum}; laplace=true)
    if laplace
        n = LinearAlgebra.checksquare(A)
        if n == 1
            return A[1, 1]
        elseif n == 2
            return A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]
        else
            return sum((-1)^(1+j) * A[1,j] * det(minor(A,j), laplace=true)
                       for j in axes(A,2))
        end
    else
        if istriu(A) || istril(A)
            return det(UpperTriangular(A))
        end
        return det(lu(A; check = false))
    end
end

function LinearAlgebra.norm(x::AbstractArray{<:RCNum}, p::Real=2)
    p = value(p)
    issym = p isa Symbolic
    if !issym && p == 2
        sqrt(sum(x->abs2(x), x))
    elseif !issym && isone(p)
        sum(abs, x)
    elseif !issym && isinf(p)
        mapreduce(abs, max, x)
    else
        sum(x->abs(x)^p, x)^inv(p)
    end
end

"""
    (a, b, islinear) = linear_expansion(t, x)

When `islinear`, return `a` and `b` such that `a * x + b == t`.
"""
function linear_expansion(t, x)
    a, b, islinear = _linear_expansion(t, x)
    x isa Num ? (wrap(a), wrap(b), islinear) : (a, b, islinear)
end
function linear_expansion(ts::AbstractArray, xs::AbstractArray)
    A = Matrix{Num}(undef, length(ts), length(xs))
    bvec = Vector{Num}(undef, length(ts))
    islinear = true
    for (i, t) in enumerate(ts)
        b = t isa Equation ? t.rhs - t.lhs : t
        for (j, x) in enumerate(xs)
            a, b, islinear = _linear_expansion(b, x)
            islinear &= islinear
            islinear || @goto FINISH
            A[i, j] = a
        end
        bvec[i] = b
    end
    @label FINISH
    return A, bvec, islinear
end
# _linear_expansion always returns `Symbolic`
function _linear_expansion(t::Equation, x)
    a₂, b₂, islinear = linear_expansion(t.rhs, x)
    islinear || return (a₂, b₂, false)
    a₁, b₁, islinear = linear_expansion(t.lhs, x)
    # t.rhs - t.lhs = 0
    return (a₂ - a₁, b₂ - b₁, islinear)
end
trival_linear_expansion(t, x) = isequal(t, x) ? (1, 0, true) : (0, t, true)

is_expansion_leaf(t) = !istree(t) || (operation(t) isa Differential)
@noinline expansion_check(op) = op isa Differential && error("The operation is a Differential. This should never happen.")
function _linear_expansion(t, x)
    t = value(t)
    t isa Symbolic || return (0, t, true)
    x = value(x)
    is_expansion_leaf(t) && return trival_linear_expansion(t, x)
    isequal(t, x) && return (1, 0, true)

    op, args = operation(t), arguments(t)
    expansion_check(op)

    if op === (+)
        a₁ = b₁ = 0
        islinear = true
        # (a₁ x + b₁) + (a₂ x + b₂) = (a₁ + a₂) x + (b₁ + b₂)
        for (i, arg) in enumerate(args)
            a₂, b₂, islinear = linear_expansion(arg, x)
            islinear || return (a₁, b₁, false)
            a₁ += a₂
            b₁ += b₂
        end
        return (a₁, b₁, true)
    elseif op === (-)
        @assert length(args) == 1
        a, b, islinear = linear_expansion(args[1], x)
        return (-a, -b, islinear)
    elseif op === (*)
        # (a₁ x + b₁) (a₂ x + b₂) = a₁ a₂ x² + (a₁ b₂ + a₂ b₁) x + b₁ b₂
        a₁ = 0
        b₁ = 1
        islinear = true
        for (i, arg) in enumerate(args)
            a₂, b₂, islinear = linear_expansion(arg, x)
            (islinear && _iszero(a₁ * a₂)) || return (a₁, b₁, false)
            a₁ = a₁ * b₂ + a₂ * b₁
            b₁ *= b₂
        end
        return (a₁, b₁, true)
    elseif op === (^)
        # (a₁ x + b₁)^(a₂ x + b₂) is linear => a₂ = 0 && (b₂ == 1 || a₁ == 0)
        a₂, b₂, islinear = linear_expansion(args[2], x)
        (islinear && _iszero(a₂)) || return (0, 0, false)
        a₁, b₁, islinear = linear_expansion(args[1], x)
        _isone(b₂) && return (a₁, b₁, islinear)
        (islinear && _iszero(a₁)) || return (0, 0, false)
        return (0, b₁^b₂, islinear)
    elseif op === (/)
        # (a₁ x + b₁)/(a₂ x + b₂) is linear => a₂ = 0
        a₂, b₂, islinear = linear_expansion(args[2], x)
        (islinear && _iszero(a₂)) || return (0, 0, false)
        a₁, b₁, islinear = linear_expansion(args[1], x)
        # (a₁ x + b₁)/b₂
        return islinear ? (a₁ / b₂, b₁ / b₂, islinear) : (0, 0, false)
    else
        for (i, arg) in enumerate(args)
            a, b, islinear = linear_expansion(arg, x)
            (_iszero(a) && islinear) || return (0, 0, false)
        end
        return (0, t, true)
    end
end

###
### Utilities
###

# Pretty much just copy-pasted from stdlib
SparseArrays.SparseMatrixCSC{Tv,Ti}(M::StridedMatrix) where {Tv<:RCNum,Ti} = _sparse(Tv,  Ti, M)
function _sparse(::Type{Tv}, ::Type{Ti}, M) where {Tv, Ti}
    nz = count(!_iszero, M)
    colptr = zeros(Ti, size(M, 2) + 1)
    nzval = Vector{Tv}(undef, nz)
    rowval = Vector{Ti}(undef, nz)
    colptr[1] = 1
    cnt = 1
    @inbounds for j = 1:size(M, 2)
        for i = 1:size(M, 1)
            v = M[i, j]
            if !_iszero(v)
                rowval[cnt] = i
                nzval[cnt] = v
                cnt += 1
            end
        end
        colptr[j+1] = cnt
    end
    return SparseMatrixCSC(size(M, 1), size(M, 2), colptr, rowval, nzval)
end
