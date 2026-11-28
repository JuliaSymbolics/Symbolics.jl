"""
    transfer_function(A, B, C, D, s; simplify = false)

Return the transfer function matrix `G(s) = C (sI - A)⁻¹ B + D` of the state-space system
`x' = Ax + Bu`, `y = Cx + Du` as a `Matrix{Num}` of size `(ny, nu)`.

`A` is an `(nx, nx)` matrix, `B` an `(nx, nu)` matrix, or a vector of length `nx` for a
single input, and `C` an `(ny, nx)` matrix, or a vector of length `nx` for a single
output. `D` is an `(ny, nu)` matrix, or a number if `ny == nu == 1`. Other sizes throw a
`DimensionMismatch`. A system without states (`A` of size `(0, 0)`) returns `D`.

Each entry is a fraction of two expanded polynomials, after clearing the denominators of
the entries of `A`, `B`, `C` and `D`. Common factors of numerator and denominator, such as
uncontrollable or unobservable modes, are only cancelled with `simplify = true`, which
calls `simplify_fractions` on each entry and can be slow. The determinants are expanded
without division, so the cost and the size of the result can grow exponentially with `nx`
for dense symbolic `A`. As in `expand`, integer coefficients use machine integers and
overflow silently, so give large numeric values as floating-point numbers.

Only the explicit form above is supported, not the descriptor form `Ex' = Ax + Bu`. If the
matrices come from symbolic differentiation (as with ModelingToolkit's
`linearize_symbolic`), the model must consist of operations that can be differentiated
symbolically.

# Examples

```jldoctest
julia> using Symbolics

julia> @variables s m c k;

julia> transfer_function([0 1; -k/m -c/m], [0, 1/m], [1, 0], 0, s)
1×1 Matrix{Num}:
 1 / (k + c*s + m*(s^2))
```

See also [`frequency_response`](@ref).
"""
function transfer_function(
        A::AbstractMatrix{<:Real}, B::AbstractVecOrMat{<:Real},
        C::AbstractVecOrMat{<:Real}, D::Union{Real, AbstractVecOrMat{<:Real}}, s::Num;
        simplify::Bool = false
    )
    nx = size(A, 1)
    Cm = C isa AbstractVector ? permutedims(collect(C)) : collect(C)
    Dm = D isa Real ? fill(D, 1, 1) : collect(D)
    ny, nu = size(Cm, 1), size(B, 2)
    size(A, 2) == nx || throw(DimensionMismatch("`A` must be square, got size $(size(A))"))
    size(B, 1) == nx ||
        throw(DimensionMismatch("`B` must have $nx rows, got size $(size(B))"))
    size(Cm, 2) == nx ||
        throw(DimensionMismatch("`C` must have $nx columns, got size $(size(C))"))
    size(Dm) == (ny, nu) ||
        throw(DimensionMismatch("`D` must have size $((ny, nu)), got $(size(D))"))

    # Row `i` of `[A B]` and of `[C D]` is multiplied by its common denominator. Then
    # `M = s*diagm(E) - A` and `G[i, j] = det([M B[:, j]; -C[i, :]' D[i, j]]) / (F[i] det(M))`.
    E, AB = _clear_denominators(Matrix{Num}([collect(A) collect(B)]))
    F, CD = _clear_denominators(Matrix{Num}([Cm Dm]))
    M = Matrix{Num}([(i == j ? E[i] * s : Num(0)) - AB[i, j] for i in 1:nx, j in 1:nx])
    detM = _expanded_det(M)
    G = Matrix{Num}(undef, ny, nu)
    for i in 1:ny, j in 1:nu
        bordered = [M AB[:, nx + j]; -permutedims(CD[i, 1:nx]) CD[i, nx + j]]
        num = _expanded_det(Matrix{Num}(bordered))
        G[i, j] = _iszero(num) ? Num(0) : num / expand(F[i] * detM)
    end
    return simplify ? map(simplify_fractions, G) : G
end

# Returns `(d, P)` with `P == d .* X` row-wise and no fractions left in `P`, so that `expand`
# can cancel like terms in the determinants.
function _clear_denominators(X::Matrix{Num})
    d = fill(Num(1), size(X, 1))
    P = similar(X)
    for i in axes(X, 1)
        fractions = map(X[i, :]) do x
            ux = unwrap(x)
            map(_fraction, isadd(ux) ? arguments(ux) : [ux])
        end
        common = Dict{SymbolicT, Int}()
        for entry in fractions, (_, den) in entry
            mergewith!(max, common, den)
        end
        for (base, pow) in common
            d[i] *= Num(base)^pow
        end
        for j in axes(X, 2)
            P[i, j] = expand(
                sum(fractions[j]) do (num, den)
                    for (base, pow) in common
                        num *= Num(base)^(pow - get(den, base, 0))
                    end
                    num
                end
            )
        end
    end
    return d, P
end

# A term as `(numerator, factors of the denominator => multiplicity)`
function _fraction(t)
    den = Dict{SymbolicT, Int}()
    isdiv(t) || return Num(t), den
    num, d = arguments(t)
    for factor in (ismul(d) ? arguments(d) : [d])
        base, pow = factor, 1
        if ispow(factor)
            b, p = arguments(factor)
            n = unwrap_const(p)
            if n isa Integer && n > 0
                base, pow = b, Int(n)
            end
        end
        _isone(base) || (den[base] = get(den, base, 0) + pow)
    end
    return Num(num), den
end

# Division-free Laplace expansion, memoised on the set of remaining columns. Zero entries
# are skipped, so sparse matrices stay cheap.
function _expanded_det(M::Matrix{Num})
    n = size(M, 1)
    memo = Dict{BitSet, Num}()
    function minor_det(row, cols)
        row > n && return Num(1)
        return get!(memo, cols) do
            total = Num[]
            for (k, j) in enumerate(cols)
                _iszero(M[row, j]) && continue
                term = M[row, j] * minor_det(row + 1, setdiff(cols, j))
                push!(total, isodd(k) ? term : -term)
            end
            expand(sum(total; init = Num(0)))
        end
    end
    return minor_det(1, BitSet(1:n))
end

"""
    frequency_response(G, s, ωs)

Evaluate the transfer function matrix `G` (see [`transfer_function`](@ref)) at `s = im*ω`
for every angular frequency `ω` in `ωs`, in `ComplexF64` arithmetic. Returns an
`Array{ComplexF64, 3}` of size `(ny, nu, length(ωs))`. The function for `G` is generated
once with [`build_function`](@ref).

`s` must be the only variable in `G`: substitute numeric values for all parameters first,
otherwise an `ArgumentError` is thrown. At a pole on the grid the entry is not finite (it
contains `Inf` or `NaN`). A common factor of numerator and denominator which vanishes on
the grid also gives `NaN`.

# Examples

```jldoctest
julia> using Symbolics

julia> @variables s T;

julia> G = transfer_function(fill(-1/T, 1, 1), [1/T], [1], 0, s);

julia> frequency_response(substitute(G, Dict(T => 2)), s, [0.5, 1.0])
1×1×2 Array{ComplexF64, 3}:
[:, :, 1] =
 0.5 - 0.5im

[:, :, 2] =
 0.2 - 0.4im
```
"""
function frequency_response(G::AbstractMatrix{Num}, s::Num, ωs::AbstractVector{<:Real})
    free = Set{SymbolicT}()
    for g in G
        union!(free, get_variables(g))
    end
    delete!(free, unwrap(s))
    isempty(free) || throw(
        ArgumentError(
            "`G` contains variables other than `$s`: $(join(free, ", ")). \
            `substitute` their values first."
        )
    )
    # NaNMath functions do not accept complex arguments
    _, evaluate! = build_function(G, s; expression = Val{false}, nanmath = false)
    response = Array{ComplexF64}(undef, size(G)..., length(ωs))
    for (k, ω) in enumerate(ωs)
        evaluate!(view(response, :, :, k), complex(0.0, ω))
    end
    return response
end
