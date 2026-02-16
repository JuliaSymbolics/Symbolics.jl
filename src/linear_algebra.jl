function nterms(t::SymbolicT)
    if iscall(t)
        return sum(nterms, arguments(t))
    else
        return 1
    end
end
nterms(t::Num) = nterms(unwrap(t))

# Soft pivoted
function sym_lu(A::AbstractMatrix{Num}; check=true)
    SINGULAR = typemax(Int)
    m, n = size(A)
    F = Matrix{Num}(undef, size(A)...)
    copyto!(F, A)
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

        if amin == SINGULAR && iszero(info)
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

function solve_for(eq::Any, var::Any; simplify=false, check=true)
    Base.depwarn("solve_for is deprecated, please use symbolic_linear_solve instead.", :solve_for)
    return symbolic_linear_solve(eq, var; simplify=simplify, check=check)
end

"""
$(TYPEDSIGNATURES)

Solve equation(s) `eqs` for a set of variables `vars`.

Assumes `length(eqs) == length(vars)`

Currently only works if all equations are linear. `check` if the expr is linear
w.r.t `vars`.

# Examples
```jldoctest
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
    x = __solve(a, b, simplify)
    simplify || return x
    if x isa SymbolicT
        return SymbolicUtils.simplify(x)
    end
    map!(SymbolicUtils.simplify, x, x)
    return x
end

function __solve(a::AbstractArray{SymbolicT}, b::AbstractArray{SymbolicT}, simplify::Bool)
    @inbounds for i in eachindex(b)
        b[i] = -b[i]
    end
    return _solve(a, b, simplify)
end
function __solve(a::SymbolicT, b::SymbolicT, simplify::Bool)
    if (a === COMMON_ZERO) || _iszero(a)
        return SymbolicUtils.Const{VartypeT}(NaN)
    end
    return a \ -b
end
__solve(a::Num, b::Num, simplify::Bool) = Num(__solve(unwrap(a), unwrap(b), simplify))

symbolic_linear_solve(eq::Equation, var::T; x...) where {T<:AbstractArray} = symbolic_linear_solve([eq], var; x...)
symbolic_linear_solve(eq::T, var::Num; x...) where {T<:AbstractArray} = first(symbolic_linear_solve(eq, [var]; x...))


function _solve(A::AbstractMatrix{SymbolicT}, b, do_simplify)
    _solve(Num.(A), b, do_simplify)
end
function _solve(A::AbstractMatrix{Num}, b::Union{AbstractArray{Num}, AbstractArray{SymbolicT}}, do_simplify)
    for i in eachindex(A)
        A[i] = SymbolicUtils.quick_cancel(A[i])
    end
    for i in eachindex(b)
        b[i] = SymbolicUtils.quick_cancel(b[i])
    end
    sol = sym_lu(A) \ b
    do_simplify ? SymbolicUtils.simplify_fractions.(sol) : sol
end

LinearAlgebra.ldiv!(A::UpperTriangular{<:Union{BasicSymbolic,RCNum}}, b::AbstractVector{<:Union{BasicSymbolic,RCNum}}, x::AbstractVector{<:Union{BasicSymbolic,RCNum}} = b) = symsub!(A, b, x)
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

LinearAlgebra.ldiv!(A::UnitLowerTriangular{<:Union{BasicSymbolic,RCNum}}, b::AbstractVector{<:Union{BasicSymbolic,RCNum}}, x::AbstractVector{<:Union{BasicSymbolic,RCNum}} = b) = symsub!(A, b, x)
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

"""
    (a, b, islinear) = linear_expansion(t, x)

When `islinear`, return `a` and `b` such that `a * x + b == t`. Instead of calling
`linear_expansion` multiple times with the same `x`, prefer using
[`Symbolics.LinearExpander`](@ref).
"""
function linear_expansion(t, x::Num)
    a, b, islin = linear_expansion(t, unwrap(x))
    Num(a), Num(b), islin
end

@inline function linear_expansion(t, x::SymbolicT)
    return LinearExpander(x)(unwrap(t))
end

function linear_expansion(ts::AbstractArray{T}, xs::AbstractArray{S}) where {T <: Union{Num, SymbolicT, Equation}, S <: Union{SymbolicT, Num}}
    ts = vec(ts)
    xs = vec(xs)
    A = Matrix{SymbolicT}(undef, length(ts), length(xs))
    bvec = Vector{SymbolicT}(undef, length(ts))
    map!(bvec, ts) do t
        if t isa Equation
            return t.rhs - t.lhs
        elseif t isa Num
            return unwrap(t)
        else
            return t
        end
    end
    for (j, x) in enumerate(xs)
        lex = LinearExpander(unwrap(x))
        for (i, t) in enumerate(bvec)
            a, resid, islin = lex(t)
            islin || return A, bvec, false
            A[i, j] = a
            bvec[i] = resid
        end
    end
    return A, bvec, true
end
@inline trivial_linear_expansion(t, x) = isequal(t, x) ? (COMMON_ONE, COMMON_ZERO, true) : (COMMON_ZERO, t, true)

is_expansion_leaf(t) = !iscall(t) || (operation(t) isa Operator)
@noinline throw_bad_expansion() = error("The operation is an Operator. This should never happen.")

const LinearExpansionResultT = Tuple{SymbolicT, SymbolicT, Bool}
"""
    $TYPEDEF

A functor which acts as the workhorse for [`Symbolics.linear_expansion`](@ref). The operation
`linear_expansion(t::Symbolics.SymbolicT, x::Symbolics.SymbolicT)` can equivalently be
phrased as `LinearExpander(x)(t)`. This functor caches data during processing, and thus
several `linear_expansion` calls with the same `x` can be rephrased as multiple calls to
`LinearExpansion(x)`. This can help avoid repeated computation and allocations.
"""
struct LinearExpander
    x::SymbolicT
    arr::Union{SymbolicT, Nothing}
    occursin_cache::Dict{SymbolicT, Bool}
    result_cache::Dict{SymbolicT, LinearExpansionResultT}
    strict::Bool
end

"""
    $TYPEDSIGNATURES

Create a `LinearExpander`. If `strict == true`, `ex` occurring in the condition of `ifelse`
will not be considered linear.
"""
function LinearExpander(ex::SymbolicT; strict = false)
    @match ex begin
        BSImpl.Term(; f, args) && if f === getindex end => begin
            LinearExpander(ex, args[1], Dict{SymbolicT, Bool}(ex => true, args[1] => true), Dict{SymbolicT, LinearExpansionResultT}(), strict)
        end
        _ => LinearExpander(ex, nothing, Dict{SymbolicT, Bool}(ex => true), Dict{SymbolicT, LinearExpansionResultT}(), strict)
    end
end
LinearExpander(ex::Num; kw...) = LinearExpander(unwrap(ex); kw...)

@inline function _linear_expansion_predicate(lex::LinearExpander, t::SymbolicT)
    (; x, arr) = lex
    isequal(t, x) && return true
    arr isa SymbolicT && isequal(t, arr) && return true
end

function _linear_expansion_occursin(lex::LinearExpander, t::SymbolicT)
    (; x, arr, occursin_cache) = lex
    cached = get(occursin_cache, t, nothing)
    if cached isa Bool
        return cached
    end
    _linear_expansion_predicate(lex, t) && return true
    iscall(t) || return false
    return occursin_cache[t] = any(Base.Fix1(_linear_expansion_occursin, lex), arguments(t))
end

function _linear_expansion_recurse(lex::LinearExpander, t::SymbolicT)
    result = get(lex.result_cache, t, nothing)
    if result isa LinearExpansionResultT
        return result
    end
    result = lex(t)
    lex.result_cache[t] = result
    return result
end

# _linear_expansion always returns `SymbolicT`
(lex::LinearExpander)(t) = (COMMON_ZERO, SConst(t), true)

function (lex::LinearExpander)(t::Equation)
    a₂, b₂, islinear = lex(t.rhs)
    islinear || return (COMMON_ZERO, COMMON_ZERO, false)
    a₁, b₁, islinear = lex(t.lhs)
    # t.rhs - t.lhs = 0
    return (a₂ - a₁, b₂ - b₁, islinear)
end

(lex::LinearExpander)(t::Num) = lex(unwrap(t))

function (lex::LinearExpander)(t::SymbolicT)
    (; x) = lex
    is_expansion_leaf(t) && return trivial_linear_expansion(t, x)
    isequal(t, x) && return (COMMON_ONE, COMMON_ZERO, true)
    _linear_expansion_occursin(lex, t) || return (COMMON_ZERO, t, true)

    @match t begin
        BSImpl.Term(; f) && if f isa Operator end => throw_bad_expansion()
        BSImpl.AddMul(; coeff, dict, variant, type, shape) => begin
            cf = Const{VartypeT}(coeff)
            if variant === SymbolicUtils.AddMulVariant.ADD
                a_buffer = SArgsT()
                b_buffer = SArgsT()
                for (k, v) in dict
                    a, b, islin = _linear_expansion_recurse(lex, k)
                    islin || return (COMMON_ZERO, COMMON_ZERO, false)
                    push!(a_buffer, a * v)
                    push!(b_buffer, b * v)
                end
                if !_iszero(coeff)
                    push!(b_buffer, cf)
                end
                a = SymbolicUtils.add_worker(VartypeT, a_buffer)
                b = SymbolicUtils.add_worker(VartypeT, b_buffer)
                return a, b, true
            else
                a = COMMON_ZERO
                b = COMMON_ZERO
                k = COMMON_ZERO
                newdict = dict
                dirty = false
                extras = SArgsT()
                for (_k, _v) in dict
                    _a, _b, islin = _linear_expansion_recurse(lex, _k)
                    islin || return (COMMON_ZERO, COMMON_ZERO, false)
                    _a_zero = _iszero(_a)
                    if _a_zero
                        isequal(_b, _k) && continue
                        if !dirty
                            newdict = copy(dict)
                            dirty = true
                        end
                        delete!(newdict, _k)
                        if SymbolicUtils.ismul(_b)
                            push!(extras, _b ^ _v)
                        else
                            newdict[_b] = _v
                        end
                        continue
                    end
                    _a_zero && continue
                    a === COMMON_ZERO || return (COMMON_ZERO, COMMON_ZERO, false)
                    _isone(_v) || return (COMMON_ZERO, COMMON_ZERO, false)
                    a = _a
                    b = _b
                    k = _k
                end
                if a === COMMON_ZERO
                    dirty || return (COMMON_ZERO, t, true)
                    b = SymbolicUtils.Mul{VartypeT}(coeff, newdict; type, shape)
                    if !isempty(extras)
                        push!(extras, b)
                        b = SymbolicUtils.mul_worker(VartypeT, extras)
                    end
                    return (COMMON_ZERO, b, true)
                end
                if !dirty
                    newdict = copy(dict)
                end
                delete!(newdict, k)
                tmp = Symbolics.Mul{VartypeT}(coeff, newdict; type, shape)
                if !isempty(extras)
                    push!(extras, tmp)
                    tmp = SymbolicUtils.mul_worker(VartypeT, extras)
                end
                return (a * tmp, b * tmp, true)
            end
        end
        BSImpl.Term(; f, args) => begin
            if f === (-)
                @assert length(args) == 1
                a, b, islin = _linear_expansion_recurse(lex, args[1])
                islin || return (COMMON_ZERO, COMMON_ZERO, false)
                return -a, -b, islin
            elseif f === (^)
                # If we are here, then both the base and exponent are non-trivial yet the
                # expression contains `x`. Thus `t` cannot be linear in `x`
                return (COMMON_ZERO, COMMON_ZERO, false)
            elseif f === getindex
                # This is a `getindex` expression that somehow contains `x`. It is not
                # identical to `x`, otherwise we would have exited earlier.
                #
                # Case 1: `x` is an indexed array and `t` is another element of the same
                # array
                arrt = args[1]
                @match x begin
                    BSImpl.Term(; f = fx, args = argsx) && if fx === getindex end => begin
                        isequal(arrt, argsx[1]) && return (COMMON_ZERO, t, true)
                    end
                    _ => nothing
                end

                indexed_t = scalarize(t)
                # when indexing a registered function/callable symbolic
                # scalarizing and indexing leads to the same symbolic variable
                # which causes a StackOverflowError without this
                isequal(t, indexed_t) && return (COMMON_ZERO, COMMON_ZERO, false)
                return _linear_expansion_recurse(lex, indexed_t)
            elseif f === ifelse
                cond, iftrue, iffalse = args
                if lex.strict && _linear_expansion_occursin(lex, cond)
                    return (COMMON_ZERO, COMMON_ZERO, false)
                end
                truea, trueb, istruelinear = _linear_expansion_recurse(lex, iftrue)
                istruelinear || return (COMMON_ZERO, t, false)
                falsea, falseb, isfalselinear = _linear_expansion_recurse(lex, iffalse)
                isfalselinear || return (COMMON_ZERO, t, false)
                a = isequal(truea, falsea) ? truea : ifelse(cond, truea, falsea)
                b = isequal(trueb, falseb) ? trueb : ifelse(cond, trueb, falseb)
                return (a, b, true)
            elseif f === LinearAlgebra.dot
                a, b = args
                sh = SymbolicUtils.shape(a)
                if !SymbolicUtils.is_array_shape(sh)
                    return _linear_expansion_recurse(lex, a * b)
                end
                if sh isa SymbolicUtils.Unknown
                    return (COMMON_ZERO, t, false)
                end
                sh = sh::SymbolicUtils.ShapeVecT
                add_buffer = SymbolicUtils.ArgsT{VartypeT}()
                sizehint!(add_buffer, length(a))
                for i in SymbolicUtils.stable_eachindex(a)
                    push!(add_buffer, LinearAlgebra.dot(a[i], b[i]))
                end
                res = SymbolicUtils.add_worker(VartypeT, add_buffer)
                return _linear_expansion_recurse(lex, res)
            else
                for (i, arg) in enumerate(args)
                    _linear_expansion_predicate(lex, arg) && return (COMMON_ZERO, COMMON_ZERO, false)
                    a, b, islin = _linear_expansion_recurse(lex, arg)
                    (_iszero(a) && islin) || return (COMMON_ZERO, COMMON_ZERO, false)
                end
                return (COMMON_ZERO, t, true)
            end
        end
        BSImpl.Div(; num, den, type, shape) => begin
            if _linear_expansion_occursin(lex, den)
                return (COMMON_ZERO, COMMON_ZERO, false)
            end
            a, b, islin = _linear_expansion_recurse(lex, num)
            islin || return (COMMON_ZERO, COMMON_ZERO, false)
            a = SymbolicUtils.Div{VartypeT}(a, den, false; type, shape)
            b = SymbolicUtils.Div{VartypeT}(b, den, false; type, shape)
            return a, b, true
        end
        BSImpl.ArrayOp(; output_idx) => begin
            if isempty(output_idx)
                # is scalar
                scal = SymbolicUtils.scalarize(t, Val{true}())::SymbolicT
                return _linear_expansion_recurse(lex, scal)
            else
                # TODO: call `_linear_expansion(t.expr, x)` and parse the result
                return COMMON_ZERO, COMMON_ZERO, false
            end
        end
    end
end

###
### Utilities
###

# Pretty much just copy-pasted from stdlib
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
