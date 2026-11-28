"""
    transfer_function(A, B, C, D, s; simplify = false)

Return the transfer function matrix `G(s) = C (sI - A)⁻¹ B + D` of the state-space system
`x' = Ax + Bu`, `y = Cx + Du` as a `Matrix{Num}` of size `(ny, nu)`.

`A` is an `(nx, nx)` matrix, `B` an `(nx, nu)` matrix, or a vector of length `nx` for a
single input, and `C` an `(ny, nx)` matrix, or a vector of length `nx` for a single
output. `D` is an `(ny, nu)` matrix, or a number if `ny == nu == 1`. Other sizes throw a
`DimensionMismatch`. A system without states (`A` of size `(0, 0)`) returns `D`.

`s` may appear in `B` and `D`, as in `B = B1 + s*B2` and `D = D1 + s*D2` for a system
driven by `u` and its derivative `u'`. The formula is applied as written, also if `A` or
`C` contain `s`; an `ArgumentError` is thrown if `sI - A` is then identically singular.
For a linearized nonlinear model, the entries contain states and inputs, which must be
substituted by their values at the operating point before calling
[`frequency_response`](@ref).

Each entry is a fraction of two expanded polynomials, after clearing the denominators of
the entries of `A`, `B`, `C` and `D`. Common factors of numerator and denominator, such as
uncontrollable or unobservable modes, are only cancelled with `simplify = true`, which
calls `simplify_fractions` on each entry and can be slow. The determinants are expanded
without division, so the cost and the size of the result can grow exponentially with `nx`
for dense symbolic `A`.

Arithmetic on constants is exact, with floating-point constants taken as the rationals they
represent. Without floating-point constants, coefficients which do not fit in an `Int` are
returned as `BigInt`s. With them, all coefficients are rounded to `Float64` at the end (to
`BigFloat` if there is a `BigFloat` constant). An `ArgumentError` is thrown for `Inf` and
`NaN` constants, and if a coefficient overflows or underflows `Float64`. Constants inside
other functions, such as `sin`, are left as given.

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

    # All arithmetic is done on polynomials with `BigInt`-based coefficients, because
    # `expand` computes in machine integers, which wrap around silently.
    atoms = Bijections.Bijection{SymbolicUtils.PolyVarT, SymbolicT}()
    floattype = Ref{Type}(Union{})
    polynomial(x) = _polynomial!(atoms, floattype, unwrap(x))

    # Row `i` of `[A B]` and of `[C D]` is multiplied by its common denominator. Then
    # `M = s*diagm(E) - A` and `G[i, j] = det([M B[:, j]; -C[i, :]' D[i, j]]) / (F[i] det(M))`.
    # converted before concatenation, which would promote rationals to floats
    E, AB = _clear_denominators(polynomial, Matrix{Num}([Num.(collect(A)) Num.(collect(B))]))
    F, CD = _clear_denominators(polynomial, Matrix{Num}([Num.(Cm) Num.(Dm)]))
    sp, o = polynomial(s), polynomial(Num(0))
    M = DP.Polynomial[(i == j ? E[i] * sp : o) - AB[i, j] for i in 1:nx, j in 1:nx]
    detM = _determinant(M, polynomial(Num(1)))
    iszero(detM) && throw(ArgumentError("`sI - A` is singular for every `$s`"))
    narrow(p) = Num(
        SymbolicUtils.from_poly(
            atoms,
            DP.Polynomial(
                SymbolicUtils.PolyCoeffT[_narrow(c, floattype[]) for c in MP.coefficients(p)],
                MP.monomials(p)
            )
        )
    )
    G = Matrix{Num}(undef, ny, nu)
    for i in 1:ny, j in 1:nu
        bordered = DP.Polynomial[M AB[:, nx + j]; -permutedims(CD[i, 1:nx]) CD[i, nx + j]]
        num = _determinant(bordered, polynomial(Num(1)))
        den = F[i] * detM
        # dividing a rational number by one would give a float
        G[i, j] = iszero(num) ? Num(0) : isone(den) ? narrow(num) : narrow(num) / narrow(den)
    end
    return simplify ? map(simplify_fractions, G) : G
end

# `x` as a polynomial with exact coefficients, where `floattype` records the type to round
# them to. Anything but `+`, `*` and non-negative integer powers becomes a variable of its
# own, recorded in `atoms`.
function _polynomial!(atoms, floattype, x::SymbolicT)
    if SymbolicUtils.isconst(x)
        v = unwrap_const(x)
        if v isa AbstractFloat
            isfinite(v) || throw(ArgumentError("the constant $v is not finite"))
            floattype[] = promote_type(floattype[], Float64, typeof(v))
        end
        c = v isa Integer ? BigInt(v) :
            v isa Union{Rational, AbstractFloat} ? Rational{BigInt}(v) : nothing
        c === nothing || return MP.polynomial(c * one(SymbolicUtils.PolyVarT))
    elseif iscall(x)
        op, args = operation(x), arguments(x)
        op === (+) && return sum(a -> _polynomial!(atoms, floattype, a), args)
        op === (*) && return prod(a -> _polynomial!(atoms, floattype, a), args)
        if op === (^)
            n = unwrap_const(args[2])
            n isa Integer && n >= 0 && return _polynomial!(atoms, floattype, args[1])^Int(n)
        end
    end
    var = SymbolicUtils.basicsymbolic_to_polyvar(Bijections.active_inv(atoms), x)
    return MP.polynomial(big(1) * var)
end

# Coefficients are rounded to `F`, or converted back to machine numbers if they fit
function _narrow(v, F)
    if F === Union{}
        q = isone(denominator(v)) ? numerator(v) : v
        fits = abs(numerator(q)) <= typemax(Int) && denominator(q) <= typemax(Int)
        return fits ? (q isa BigInt ? Int(q) : Rational{Int}(q)) : q
    end
    x = F(v)
    isfinite(x) && iszero(x) == iszero(v) || throw(
        ArgumentError("the coefficient $(BigFloat(v)) cannot be represented as a `$F`")
    )
    return x
end

# Returns `(d, P)` with `P == d .* X` row-wise and no fractions left in `P`, so that like
# terms cancel in the determinants.
function _clear_denominators(polynomial, X::Matrix{Num})
    d = DP.Polynomial[polynomial(Num(1)) for _ in axes(X, 1)]
    P = Matrix{DP.Polynomial}(undef, size(X))
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
            d[i] *= polynomial(base)^pow
        end
        for j in axes(X, 2)
            P[i, j] = sum(fractions[j]) do (num, den)
                term = polynomial(num)
                for (base, pow) in common
                    term *= polynomial(base)^(pow - get(den, base, 0))
                end
                term
            end
        end
    end
    return d, P
end

# A term as `(numerator, factors of the denominator => multiplicity)`
function _fraction(t)
    den = Dict{SymbolicT, Int}()
    isdiv(t) || return t, den
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
    return num, den
end

# Division-free Laplace expansion, memoised on the set of remaining columns. Zero entries
# are skipped, so sparse matrices stay cheap.
function _determinant(M::Matrix{DP.Polynomial}, unit)
    n = size(M, 1)
    memo = Dict{BitSet, DP.Polynomial}()
    function minor_det(row, cols)
        row > n && return unit
        return get!(memo, cols) do
            total = zero(unit)
            for (k, j) in enumerate(cols)
                iszero(M[row, j]) && continue
                term = M[row, j] * minor_det(row + 1, setdiff(cols, j))
                total += isodd(k) ? term : -term
            end
            total
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
otherwise an `ArgumentError` is thrown. `substitute` computes with machine integers, which
wrap around silently, so use floating-point numbers for large parameter values. At a pole on the grid the entry is not finite (it
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
