function nterms(t)
    if iscall(t)
        return reduce(+, map(nterms, arguments(t)), init=0)
    else
        return 1
    end
end

# Soft pivoted
# Note: we call this function with a matrix of Union{SymbolicUtils.Symbolic, Any}
function sym_lu(A; check=true)
    SINGULAR = typemax(Int)
    m, n = size(A)
    F = map(x->x isa RCNum ? x : Num(x), A)
    minmn = min(m, n)
    p = Vector{LinearAlgebra.BlasInt}(undef, minmn)
    info = 0
    for k = 1:minmn
        kp = k
        amin = SINGULAR
        for i in k:m
            absi = _iszero(F[i, k]) ? SINGULAR : nterms(F[i,k])
            if absi < amin
                kp = i
                amin = absi
            end
        end

        p[k] = kp

        if amin == SINGULAR && !(amin isa Symbolic) && (amin isa Number) && iszero(info)
            info = k
        end

        # swap
        for j in 1:n
            F[k, j], F[kp, j] = F[kp, j], F[k, j]
        end

        for i in k+1:m
            F[i, k] = F[i, k] / F[k, k]
        end
        for j = k+1:n
            for i in k+1:m
                F[i, j] = F[i, j] - F[i, k] * F[k, j]
            end
        end
    end
    check && LinearAlgebra.checknonsingular(info)
    LU(F, p, convert(LinearAlgebra.BlasInt, info))
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

function solve_for(eq::Any, var::Any; simplify=false, check=true)
    Base.depwarn("solve_for is deprecated, please use symbolic_linear_solve instead.", :solve_for, force=true)
    return symbolic_linear_solve(eq, var; simplify=simplify, check=check)
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

julia> Symbolics.symbolic_linear_solve(x + y ~ 0, x)
-y

julia> Symbolics.symbolic_linear_solve([x + y ~ 0, x - y ~ 2], [x, y])
2-element Vector{Float64}:
  1.0
 -1.0
```
"""
function symbolic_linear_solve(eq, var; simplify=false, check=true) # scalar case
    # simplify defaults for `false` as canonicalization should handle most of
    # the cases.
    a, b, islinear = linear_expansion(eq, var)
    check && @assert islinear
    islinear || return nothing
    # a * x + b = 0
    if eq isa AbstractArray && var isa AbstractArray
        x = _solve(a, -b, simplify)
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
symbolic_linear_solve(eq::Equation, var::T; x...) where {T<:AbstractArray} = symbolic_linear_solve([eq], var; x...)
symbolic_linear_solve(eq::T, var::Num; x...) where {T<:AbstractArray} = first(symbolic_linear_solve(eq, [var]; x...))


function _solve(A::AbstractMatrix, b::AbstractArray, do_simplify)
    A = Num.(SymbolicUtils.quick_cancel.(A))
    b = Num.(SymbolicUtils.quick_cancel.(b))
    sol = value.(sym_lu(A) \ b)
    do_simplify ? SymbolicUtils.simplify_fractions.(sol) : sol
end

LinearAlgebra.ldiv!(A::UpperTriangular{<:Union{Symbolic,RCNum}}, b::AbstractVector{<:Union{Symbolic,RCNum}}, x::AbstractVector{<:Union{Symbolic,RCNum}} = b) = symsub!(A, b, x)
function symsub!(A::UpperTriangular, b::AbstractVector, x::AbstractVector = b)
    LinearAlgebra.require_one_based_indexing(A, b, x)
    n = size(A, 2)
    if !(n == length(b) == length(x))
        throw(DimensionMismatch("second dimension of left hand side A, $n, length of output x, $(length(x)), and length of right hand side b, $(length(b)), must be equal"))
    end
    @inbounds for j in n:-1:1
        _iszero(A.data[j,j]) && throw(SingularException(j))
        xj = x[j] = b[j] / A.data[j,j]
        for i in j-1:-1:1
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
    @inbounds for j in 1:n
        xj = x[j] = b[j]
        for i in j+1:n
            sub = _isone(xj) ? A.data[i,j] : A.data[i,j] * xj
            if !_iszero(sub)
                b[i] -= sub
            end
        end
    end
    x
end

minor(B, j) = @view B[2:end, 1:size(B,2) .!= j]
minor(B, i, j) = @view B[1:size(B,1) .!= i, 1:size(B,2) .!= j]
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

LinearAlgebra.inv(A::AbstractMatrix{<:RCNum}; laplace=true) = _invl(A; laplace=laplace)
LinearAlgebra.inv(A::StridedMatrix{<:RCNum}; laplace=true) = _invl(A; laplace=laplace)

function _invl(A::AbstractMatrix{<:RCNum}; laplace=true)
    if laplace
		n = LinearAlgebra.checksquare(A)
        A⁻¹ = similar(A)
        idet = 1/det(A; laplace=true)
	if n==1
	    A⁻¹[1,1] = idet
	    return A⁻¹
	end
        for i=1:size(A,1)
            for j = 1:size(A,1)
                A⁻¹[i,j] = (-1)^(i+j)*det(minor(A, j, i); laplace=true)*idet
            end
        end
        return A⁻¹
    else
        if istriu(A) || istril(A)
            return inv(UpperTriangular(A))
        end
        return inv(lu(A; check = false))
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
trivial_linear_expansion(t, x) = isequal(t, x) ? (1, 0, true) : (0, t, true)

is_expansion_leaf(t) = !iscall(t) || (operation(t) isa Operator)
@noinline expansion_check(op) = op isa Operator && error("The operation is an Operator. This should never happen.")
function _linear_expansion(t, x)
    t = value(t)
    t isa Symbolic || return (0, t, true)
    x = value(x)
    is_expansion_leaf(t) && return trivial_linear_expansion(t, x)
    isequal(t, x) && return (1, 0, true)

    op, args = operation(t), arguments(t)
    expansion_check(op)

    if iscall(x) && operation(x) == getindex
        arrx, idxsx... = arguments(x)
    else
        arrx = nothing
        idxsx = nothing
    end

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
    elseif op === getindex
        arrt, idxst... = arguments(t)
        isequal(arrt, arrx) && return (0, t, true)

        indexed_t = Symbolics.scalarize(arrt)[idxst...]
        # when indexing a registered function/callable symbolic
        # scalarizing and indexing leads to the same symbolic variable
        # which causes a StackOverflowError without this
        isequal(t, indexed_t) && return (0, t, true)
        return linear_expansion(Symbolics.scalarize(arrt)[idxst...], x)
    else
        for (i, arg) in enumerate(args)
            isequal(arg, arrx) && return (0, 0, false)
            if symbolic_type(arg) == NotSymbolic()
                arg isa AbstractArray || continue
                _occursin_array(x, arrx, arg) && return (0, 0, false)
                continue
            end
            a, b, islinear = linear_expansion(arg, x)
            (_iszero(a) && islinear) || return (0, 0, false)
        end
        return (0, t, true)
    end
end

"""
    _occursin_array(sym, arrsym, arr)

Check if `sym` (or, if `sym` is an element of an array symbolic, the array symbolic
`arrsym`) occursin in the non-symbolic array `arr`.
"""
function _occursin_array(sym, arrsym, arr)
    for el in arr
        if symbolic_type(el) == NotSymbolic()
            return el isa AbstractArray && _occursin_array(sym, arrsym, el)
        else
            return occursin(sym, el) || occursin(arrsym, el)
        end
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
    @inbounds for j in 1:size(M, 2)
        for i in 1:size(M, 1)
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
