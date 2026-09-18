using BenchmarkTools, Symbolics, Random, LinearAlgebra
using Symbolics: linear_expansion

const SUITE = BenchmarkGroup()

# ── matrix / linear-algebra ───────────────────────────────────────────────────

vars = @variables a, b, c, d, e, f, g, h, i

X = [0 b c;
     d e f;
     g h i]

F = lu(X)
FbyX = F \ X

x = (f + ((((g*(c^2)*(e^2)) / d - e*h*(c^2)) / b + (-c*e*f*g) / d + c*e*i) /
          (i + ((c*e*g) / d - c*h) / b + (-f*g) / d) - c*e) / b +
     ((g*(f^2)) / d + ((-c*e*f*g) / d + c*f*h) / b - f*i) /
     (i + ((c*e*g) / d - c*h) / b + (-f*g) / d)) / d

o = (d + (e*((c*(g + (-d*g) / d)) / (i + (-c*(h + (-e*g) / d)) / b + (-f*g) / d))) / b +
     (-f*(g + (-d*g) / d)) / (i + (-c*(h + (-e*g) / d)) / b + (-f*g) / d)) / d

SUITE["iszero/1"] = @benchmarkable iszero($((b*(h + (-e*g) / d)) / b + (e*g) / d - h))
SUITE["isone/1"]  = @benchmarkable FbyX == I
SUITE["iszero/2"] = @benchmarkable iszero(x)
SUITE["isone/2"]  = @benchmarkable isone(o)
SUITE["lu"]       = @benchmarkable lu(X)
SUITE["_solve"]   = @benchmarkable Symbolics._solve(X, X, true)

# ── random expression tree helpers ────────────────────────────────────────────

const _UNARY_OPS = [sin, cos, exp]

# Recursively build a random symbolic expression tree of the given depth.
# All leaf nodes are drawn from `atoms`. The RNG must be seeded before calling
# to ensure reproducibility.
function _random_expr(rng::AbstractRNG, atoms::AbstractVector, depth::Int)
    depth == 0 && return rand(rng, atoms)
    r = rand(rng)
    if r < 0.25      # unary: sin / cos / exp
        return rand(rng, _UNARY_OPS)(_random_expr(rng, atoms, depth - 1))
    elseif r < 0.50  # n-ary sum
        nargs = rand(rng, 2:4)
        return +([_random_expr(rng, atoms, depth - 1) for _ in 1:nargs]...)
    elseif r < 0.75  # n-ary product
        nargs = rand(rng, 2:3)
        return *([_random_expr(rng, atoms, depth - 1) for _ in 1:nargs]...)
    else             # binary difference
        return _random_expr(rng, atoms, depth - 1) - _random_expr(rng, atoms, depth - 1)
    end
end

# Build a linear expression in `var` with `depth` nested levels of multiplication:
#   coeff_d * (... (coeff_1 * var + b_1) ...) + b_d
# Symbolics distributes the products, yielding an Add with d+1 terms where the linear Mul
# has d+1 distinct complex factors. linear_expansion must iterate over every factor to
# locate var, so runtime scales with depth.
function _random_linear_expr_deep(
        rng::AbstractRNG, coeff_atoms::AbstractVector, var, depth::Int)
    depth == 0 && return var
    coeff      = _random_expr(rng, coeff_atoms, 2)
    const_term = _random_expr(rng, coeff_atoms, 2)
    return coeff * _random_linear_expr_deep(rng, coeff_atoms, var, depth - 1) + const_term
end

# Build n expressions (one per unknown), each linear in all variables in `unknowns`.
function _random_linear_system(
        rng::AbstractRNG, coeff_atoms::AbstractVector, unknowns::AbstractVector, depth::Int)
    n = length(unknowns)
    map(1:n) do _
        constant    = _random_expr(rng, coeff_atoms, depth)
        linear_part = sum(_random_expr(rng, coeff_atoms, depth) * u for u in unknowns)
        constant + linear_part
    end
end

# Build a linear expression in `var` with `depth` nested ifelse layers:
#   ifelse(cond_d, tc_d * inner + b_d, fc_d * inner + d_d)  where inner = depth-1 version.
# Unlike plain multiplication, Symbolics does NOT distribute `tc * ifelse(...)`, so the
# nested ifelse structure is preserved in the expression tree. linear_expansion recurses
# through each layer — into both branches of every ifelse — before reaching var.
function _random_ifelse_linear_expr_deep(
        rng::AbstractRNG, coeff_atoms::AbstractVector, var, depth::Int)
    depth == 0 && return var
    inner       = _random_ifelse_linear_expr_deep(rng, coeff_atoms, var, depth - 1)
    cond        = rand(rng, coeff_atoms) > rand(rng, coeff_atoms)
    true_coeff  = _random_expr(rng, coeff_atoms, 2)
    true_const  = _random_expr(rng, coeff_atoms, 2)
    false_coeff = _random_expr(rng, coeff_atoms, 2)
    false_const = _random_expr(rng, coeff_atoms, 2)
    return ifelse(cond, true_coeff * inner + true_const, false_coeff * inner + false_const)
end

# Build n expressions (one per unknown) each containing per-unknown ifelse terms linear
# in all variables in `unknowns`.
function _random_ifelse_linear_system(
        rng::AbstractRNG, coeff_atoms::AbstractVector, unknowns::AbstractVector, depth::Int)
    n = length(unknowns)
    map(1:n) do _
        constant = _random_expr(rng, coeff_atoms, depth)
        terms = map(unknowns) do u
            cond        = rand(rng, coeff_atoms) > rand(rng, coeff_atoms)
            true_coeff  = _random_expr(rng, coeff_atoms, depth)
            true_const  = _random_expr(rng, coeff_atoms, depth)
            false_coeff = _random_expr(rng, coeff_atoms, depth)
            false_const = _random_expr(rng, coeff_atoms, depth)
            ifelse(cond, true_coeff * u + true_const, false_coeff * u + false_const)
        end
        constant + sum(terms)
    end
end

# ── symbolic AD: expand_derivatives ──────────────────────────────────────────

@variables p q

const _ad_atoms = [p, q]
const _Dp = Differential(p)

const SUITE_AD = SUITE["AD"] = BenchmarkGroup()

for depth in (3, 5, 7, 9, 11)
    expr = _random_expr(Xoshiro(42), _ad_atoms, depth)
    SUITE_AD["expand_derivatives/depth=$depth"] = @benchmarkable expand_derivatives($(_Dp(expr)))
    SUITE_AD["dstar_derivative/depth=$depth"] = @benchmarkable Symbolics.dstar_derivative($expr, $p)
end

# ── symbolic AD: dstar_jacobian vs jacobian ───────────────────────────────────

@variables dx dy dz dw shx shy shz

# shared subexpression chain: p=x+y, q=x-y, r=pq, s=p/q, t=r+s, u=rs
_dj_p = dx + dy; _dj_q = dx - dy; _dj_r = _dj_p * _dj_q; _dj_s = _dj_p / _dj_q
_dj_t = _dj_r + _dj_s; _dj_u = _dj_r * _dj_s

# spherical harmonics Y(l, m) for l < max_l — heavily shared polynomial structure
function _spherical_harmonics(max_l::Integer, x, y, z)
    Pc = Dict{Tuple{Int,Int}, Any}(); Cc = Dict{Int, Any}(); Sc = Dict{Int, Any}()
    P(l, m) = get!(Pc, (l, m)) do
        if l == 0 && m == 0; 1.0
        elseif l == m; (1 - 2m) * P(m - 1, m - 1)
        elseif l == m + 1; (2m + 1) * z * P(m, m)
        else; ((2l - 1) / (l - m)) * z * P(l - 1, m) - ((l + m - 1) / (l - m)) * P(l - 2, m)
        end
    end
    C(m) = get!(Cc, m) do; m == 0 ? 1 : x * S(m - 1) + y * C(m - 1); end
    S(m) = get!(Sc, m) do; m == 0 ? 0 : x * C(m - 1) - y * S(m - 1); end
    factorial_approx(n) = sqrt(2π * n) * (n / ℯ * sqrt(n * sinh(1 / n) + 1 / (810 * n^6)))^n
    N(l, m) = m == 0 ? sqrt(2l + 1 / (4π)) : sqrt((2l + 1) / 2π * factorial_approx(l - m) / factorial_approx(l + m))
    Y(l, m) = m < 0 ? N(l, -m) * P(l, -m) * S(-m) : N(l, m) * P(l, m) * C(m)
    unwrap.([Y(l, m) for l in 0:max_l-1 for m in -l:l])
end

const SUITE_DJ = SUITE["dstar_jacobian"] = BenchmarkGroup()

# "sparse"/"mixed": roots contain subtrees with no vars — edges into them can
# never reach a requested var, exercising D*'s dead-edge handling
const _DJ_WORKLOADS = [
    ("simple", unwrap.([sin(dx)*dx^2 + dx]),                    unwrap.([dx])),
    ("shared", unwrap.([(dx^2+dy^2)^2, (dx^2+dy^2)^2*dy]),      unwrap.([dx,dy])),
    ("pq",     unwrap.([_dj_t, _dj_u, _dj_t*_dj_u]),            unwrap.([dx,dy])),
    ("sparse", unwrap.([dx*(dy^3 + sin(dy) + cos(dy*dz)) + dx*(dz^2 + exp(dz*dw)) + sin(dx)*(dy + dz),
                        dx*log(dy + dz) + dx^2*(dw + dy)]),     unwrap.([dx])),
    ("mixed",  unwrap.([dx^2 + sin(dz), dy*cos(dz), dx*dy]),    unwrap.([dx,dy])),
    ("sh4",    _spherical_harmonics(4, shx, shy, shz),          unwrap.([shx,shy,shz])),
    ("sh5",    _spherical_harmonics(5, shx, shy, shz),          unwrap.([shx,shy,shz])),
    ("sh8",    _spherical_harmonics(8, shx, shy, shz),          unwrap.([shx,shy,shz])),
]

for (name, roots, vars) in _DJ_WORKLOADS
    SUITE_DJ["$name/dstar"] = @benchmarkable Symbolics.dstar_jacobian($roots, $vars)
    SUITE_DJ["$name/jacobian"] = @benchmarkable jacobian($roots, $vars)
end

# ── linear_expansion ──────────────────────────────────────────────────────────

@variables r           # scalar unknown
@variables s[1:4]     # coefficient atoms
@variables w[1:64]    # vector unknowns (for matrix benchmarks)

const _le_unknown  = r
const _le_coeffs   = collect(s)
const _le_unknowns = collect(w)

const SUITE_LE = SUITE["linear_expansion"] = BenchmarkGroup()

# scalar: `coeff * linear_expr(var)` nested depth levels deep; each level adds one more
# factor to the linear Mul that linear_expansion must iterate through.
for depth in (2, 8, 32, 128, 512)
    expr = _random_linear_expr_deep(Xoshiro(42), _le_coeffs, _le_unknown, depth)
    SUITE_LE["scalar/depth=$depth"] = @benchmarkable linear_expansion($expr, $r)
end

# matrix: n×n system (n equations linear in n unknowns); coefficient depth fixed at 2
for n in (4, 8, 16, 32, 64)
    unknowns = _le_unknowns[1:n]
    eqs      = _random_linear_system(Xoshiro(42), _le_coeffs, unknowns, 2)
    SUITE_LE["matrix/n=$n"] = @benchmarkable linear_expansion($eqs, $unknowns)
end

# ifelse/scalar: `tc * ifelse_expr(var)` nested depth levels deep.
# Symbolics does not distribute multiplication into ifelse, so the tree depth is preserved.
# linear_expansion recurses through each ifelse layer in both branches before reaching var.
for depth in (2, 8, 32, 128, 512)
    expr = _random_ifelse_linear_expr_deep(Xoshiro(42), _le_coeffs, _le_unknown, depth)
    SUITE_LE["ifelse/scalar/depth=$depth"] = @benchmarkable linear_expansion($expr, $r)
end

# ifelse/matrix: n equations where each per-unknown term is an ifelse expression.
for n in (4, 8, 16, 32)
    unknowns = _le_unknowns[1:n]
    eqs      = _random_ifelse_linear_system(Xoshiro(42), _le_coeffs, unknowns, 2)
    SUITE_LE["ifelse/matrix/n=$n"] = @benchmarkable linear_expansion($eqs, $unknowns)
end
