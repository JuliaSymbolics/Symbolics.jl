using Symbolics
using Test
using Latexify
using LaTeXStrings
using ReferenceTests
import SymbolicUtils as SU
using Symbolics: VartypeT

using DomainSets: Interval, ClosedInterval

const LatexifyExt = Base.get_extension(Symbolics, :SymbolicsLatexifyExt)
const _toexpr = LatexifyExt._toexpr

@variables x y z u(x) dx h[1:10, 1:10] hh(x, y)[1:10, 1:10] gg(x, y)[1:10, 1:10] [latexwrapper = string]
@variables AA(x) [latexwrapper = string] X₁(x) [latexwrapper = string]
@variables a[1:10]
Dx = Differential(x)
Dy = Differential(y)

# issue 260
@test _toexpr(3 * x / y) == :((3x) / y)
@test _toexpr(3 * x^y) == :(3x^y)

@test_reference "latexify_refs/inverse.txt" latexify(x^-1)

@test_reference "latexify_refs/integral1.txt" latexify(Integral(y in Interval(0, 1))(x))
@test_reference "latexify_refs/integral2.txt" latexify(Integral(y in Interval(-Inf, Inf))(u^2))
@test_reference "latexify_refs/integral3.txt" latexify(Integral(y in Interval(-z, u))(x^2))

@test_reference "latexify_refs/frac1.txt" latexify((z + x * y^-1) / sin(z))
@test_reference "latexify_refs/frac2.txt"  latexify((3x - 7y * z^23) * (z - z^2) / x)

@test_reference "latexify_refs/minus1.txt" latexify(x - y)
@test_reference "latexify_refs/minus2.txt" latexify(x - y * z)
@test_reference "latexify_refs/minus3.txt" latexify(sin(x + y - z))

@test_reference "latexify_refs/unary_minus1.txt" latexify(-y)
@test_reference "latexify_refs/unary_minus2.txt" latexify(-y * z)

@test_reference "latexify_refs/derivative1.txt" latexify(Dx(y))
@test_reference "latexify_refs/derivative2.txt" latexify(Dx(u))
@test_reference "latexify_refs/derivative3.txt" latexify(Dx(x^2 + y^2 + z^2))
@test_reference "latexify_refs/derivative4.txt" latexify(Dy(u))
@test_reference "latexify_refs/derivative5.txt" latexify(Dx(Dy(Dx(y))))

# issue #1979: exact strings (hand-written from the issue's requested forms).
@testset "latexify derivatives/integrals (#1979)" begin
    @variables t x(t)
    D = Differential(t)
    # Issue examples with mult_symbol override: no leak into the differential form.
    @test latexify(D(x); env = :raw, mult_symbol = "\\cdot").s ==
        "\\frac{\\mathrm{d}x\\left( t \\right)}{\\mathrm{d}t}"
    @test latexify(D(D(x)); env = :raw, mult_symbol = "\\cdot").s ==
        "\\frac{\\mathrm{d}^{2}x\\left( t \\right)}{\\mathrm{d}t^{2}}"

    @variables x y a[1:3]
    Dx = Differential(x)
    I = Integral(x in Interval(0, 1))

    # Compound integrand: differential after integrand; product stays grouped.
    @test latexify(I(x * y); env = :raw).s ==
        "\\int_{0}^{1} ~ \\left( x ~ y \\right) ~ \\mathrm{d}x"
    # Sum integrand stays parenthesised (not flattened to ∫x + y).
    @test latexify(I(x + y); env = :raw).s ==
        "\\int_{0}^{1} ~ \\left( x + y \\right) ~ \\mathrm{d}x"
    # Indexed variable keeps recipe index=:subscript; no redundant merge parens.
    @test latexify(I(a[1]); env = :raw).s ==
        "\\int_{0}^{1} ~ a_{1} ~ \\mathrm{d}x"
    # Float coefficient keeps FancyNumberFormatter(5) rounding.
    @test latexify(I(0.123456789 * x); env = :raw).s ==
        "\\int_{0}^{1} ~ \\left( 0.12346 ~ x \\right) ~ \\mathrm{d}x"
    # Power of an operator-form derivative keeps the derivative parenthesised as base.
    @test latexify((Dx(x + y))^2; env = :raw).s ==
        "\\left( \\frac{\\mathrm{d}}{\\mathrm{d}x} ~ \\left( x + y \\right) \\right)^{2}"
    # Integral upper limit closes with `}` (not `)`), and mult_symbol does not leak.
    @test latexify(I(y); env = :raw, mult_symbol = "\\cdot").s ==
        "\\int_{0}^{1} ~ y ~ \\mathrm{d}x"
end

# issue #1694: bare Integral recipe and differential form.
@testset "latexify bare Integral (#1694)" begin
    @variables x a b y
    I1 = Integral(x in ClosedInterval(a, b))
    @test latexify(I1; env = :raw).s == "\\int_{a}^{b} ~ \\mathrm{d}x"
    @test latexify(I1(x^2 + 2 * x + 1); env = :raw).s ==
        "\\int_{a}^{b} ~ \\left( 1 + 2 ~ x + x^{2} \\right) ~ \\mathrm{d}x"
    @test latexify(y ~ I1(x^2 + 2 * x + 1); env = :raw).s ==
        "y = \\int_{a}^{b} ~ \\left( 1 + 2 \\cdot x + x^{2} \\right) ~ \\mathrm{d}x"
end

@test_reference "latexify_refs/stable_mul_ordering1.txt" latexify(x * y)
@test_reference "latexify_refs/stable_mul_ordering2.txt" latexify(y * x)

@test_reference "latexify_refs/equation1.txt" latexify(x ~ y + z)
@test_reference "latexify_refs/equation2.txt" latexify(x ~ Dx(y + z))

# The relative order of the two degree-1 terms AA(x) and X₁(x) in this sum is
# decided by SymbolicUtils' hash-based tie-break, which differs across Julia
# versions, so no single reference string passes on all of them. Compare the set
# of additive terms instead of the exact rendered string. (was equation5.txt)
let
    body = strip(replace(string(latexify(AA^2 + AA + 1 + X₁)),
                         "\\begin{equation}" => "", "\\end{equation}" => ""))
    @test sort(strip.(split(body, " + "))) == sort([
        "1",
        "AA\\left( x \\right)",
        "X_1\\left( x \\right)",
        "\\left( AA\\left( x \\right) \\right)^{2}",
    ])
end

@test_reference "latexify_refs/equation_vec1.txt" latexify(
    [
        x ~ y + z
        y ~ x - 3z
    ]
)
@test_reference "latexify_refs/equation_vec2.txt" latexify(
    [
        Dx(u) ~ z
        Dx(y) ~ y * x
    ]
)

@variables s p(s)[1:2] q(s)[1:2] A[1:2, 1:2]
@test_reference "latexify_refs/equation_vec_array.txt" latexify([q ~ A * p])

@test_reference "latexify_refs/complex1.txt" latexify(x^2 - y^2 + 2im * x * y)
@test_reference "latexify_refs/complex2.txt" latexify(3im * x)
@test_reference "latexify_refs/complex3.txt" latexify(1 - x + (1 + 2x) * im; imaginary_unit = "\\mathbb{i}")
@test_reference "latexify_refs/complex4.txt" latexify(im * SU.Term{VartypeT}(sqrt, [2]; type = Real, shape = []))

@syms c
@test_reference "latexify_refs/complex5.txt" latexify((3 + im / im)c)

@test_reference "latexify_refs/indices1.txt" latexify(h[10, 10])
@test_reference "latexify_refs/indices2.txt" latexify(h[10, 10], index = :bracket)

# test for https://github.com/JuliaSymbolics/Symbolics.jl/issues/1167
# note these tests need updating if/when https://github.com/korsbo/Latexify.jl/issues/331 is fixed
@test_reference "latexify_refs/indices3.txt" latexify(hh[10, 10])
@test_reference "latexify_refs/indices4.txt" latexify(gg[10, 10])

@test_reference "latexify_refs/indices5.txt" latexify(a'a)

@variables f(..)
@test_reference "latexify_refs/call_with_metadata.txt" latexify(f)

# The `_toexpr_metadata`/`_toexpr_op` hooks live in `Symbolics`, so downstream code can
# extend them via `import Symbolics` without `Base.get_extension`. The Latexify-coupled
# helpers (`_toexpr_plain`, `default_latex_wrapper`) come from the loaded extension.
struct LatexHookCtx end
function Symbolics._toexpr_metadata(O, ::Type{LatexHookCtx}, val; latexwrapper = LatexifyExt.default_latex_wrapper)
    inner = LatexifyExt._toexpr_plain(O; latexwrapper)
    return Expr(:call, :_textbf, inner)
end
expr = SU.setmetadata(x + y, LatexHookCtx, true)
@test occursin("\\textbf", latexify(expr).s)

avgf(x) = x
function Symbolics._toexpr_op(::typeof(avgf), args; latexwrapper = LatexifyExt.default_latex_wrapper)
    inner = LatexifyExt._toexpr_plain(args[1]; latexwrapper)
    inner_s = strip(latexify(inner).s, '\$')
    return LaTeXString("\\langle " * inner_s * " \\rangle")
end
avg_expr = SU.term(avgf, x + y)
@test occursin("\\langle", latexify(avg_expr).s)

@test !occursin("identity", latexify(Num(π))) # issue #1254

# issue #1820: hasmetadata should not be called on Vector arguments in getindex
@testset "getindex with vector argument (#1820)" begin
    @variables kp
    # Create an expression that indexes a literal vector of symbolic expressions.
    # This mimics what happens in piecewise/ifelse expressions with vector results.
    vec_expr = Symbolics.value.(Symbolics.Num[2kp, 3kp^2, kp + 1])
    getindex_term = SU.term(getindex, vec_expr, 1; type=SU.SymReal)
    # This should not throw a MethodError about hasmetadata on Vector
    @test_nowarn latexify(getindex_term)
end

@testset "ifelse_eager / ifelse_branching render" begin
    @variables a b
    @test_nowarn latexify(ifelse_eager(a > 0, a^2, 1 / a))
    @test_nowarn latexify(ifelse_branching(a > 0, a^2, 1 / a))
end
