using Symbolics
using Test
using Latexify
using LaTeXStrings
using ReferenceTests
import SymbolicUtils as SU
using Symbolics: VartypeT

using DomainSets: Interval

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

@test_reference "latexify_refs/integral1.txt" latexify(Integral(dx in Interval(0, 1))(x))
@test_reference "latexify_refs/integral2.txt" latexify(Integral(dx in Interval(-Inf, Inf))(u^2))
@test_reference "latexify_refs/integral3.txt" latexify(Integral(dx in Interval(-z, u))(x^2))

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

    # Compound integrand: default mult_symbol stays "~"; product stays grouped.
    @test latexify(I(x * y); env = :raw).s ==
        "\\int_{0}^{1} ~ x ~ \\left( x ~ y \\right)"
    # Sum integrand stays parenthesised (not flattened to ∫x + y).
    @test latexify(I(x + y); env = :raw).s ==
        "\\int_{0}^{1} ~ x ~ \\left( x + y \\right)"
    # Indexed variable keeps recipe index=:subscript; no redundant merge parens.
    @test latexify(I(a[1]); env = :raw).s ==
        "\\int_{0}^{1} ~ x ~ a_{1}"
    # Float coefficient keeps FancyNumberFormatter(5) rounding.
    @test latexify(I(0.123456789 * x); env = :raw).s ==
        "\\int_{0}^{1} ~ x ~ \\left( 0.12346 ~ x \\right)"
    # Power of an operator-form derivative keeps the derivative parenthesised as base.
    @test latexify((Dx(x + y))^2; env = :raw).s ==
        "\\left( \\frac{\\mathrm{d}}{\\mathrm{d}x} ~ \\left( x + y \\right) \\right)^{2}"
    # Integral upper limit closes with `}` (not `)`), and mult_symbol does not leak.
    @test latexify(I(y); env = :raw, mult_symbol = "\\cdot").s ==
        "\\int_{0}^{1} ~ x ~ y"
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

# Indexed dependent arrays (issues #1167 / #1526): subscript before call args, no escaped `\_`
@test_reference "latexify_refs/indices3.txt" latexify(hh[10, 10])
@test_reference "latexify_refs/indices4.txt" latexify(gg[10, 10])

@testset "indexed dependent array latexify (#1526)" begin
    @variables t u(t)[1:3]
    s = string(latexify(u[1]))
    @test occursin(raw"u_{1}\left( t \right)", s)
    @test !occursin(raw"\_", s)
    @test !occursin(raw"u\left( t \right)_{1}", s)
end

@testset "indexed derivative operand scope (#1526)" begin
    @variables t y(t) u(t)[1:3]
    D = Differential(t)
    @test String(latexify(D(u[1]) * y; env = :raw)) ==
        raw"\frac{\mathrm{d}u_{1}\left( t \right)}{\mathrm{d}t} ~ y\left( t \right)"
    @test String(latexify(D(D(u[1])) * y; env = :raw)) ==
        raw"\frac{\mathrm{d}^{2}u_{1}\left( t \right)}{\mathrm{d}t^{2}} ~ y\left( t \right)"
    @test String(latexify(D(u[1]) * y; env = :raw, index = :bracket)) ==
        raw"\frac{\mathrm{d}u\left[1\right]\left( t \right)}{\mathrm{d}t} ~ y\left( t \right)"
    @test String(latexify(D(D(u[1])) * y; env = :raw, index = :bracket)) ==
        raw"\frac{\mathrm{d}^{2}u\left[1\right]\left( t \right)}{\mathrm{d}t^{2}} ~ y\left( t \right)"
end

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

# issue #956: latexwrapper `_`/`^` must not be escaped when latexifying a Num
@testset "latexwrapper raw LaTeX (#956)" begin
    @variables x t
    @variables w0 [latexwrapper = s -> raw"\omega_{0}"]
    @variables vx(x, t) [latexwrapper = s -> "v_{x}"]
    ex = vx + w0^2
    eq = vx ~ w0

    for s in (string(latexify(ex)), repr(MIME"text/latex"(), ex))
        @test occursin(raw"v_{x}\left( x, t \right) + \omega_{0}^{2}", s)
        @test !occursin(raw"\_{", s)
    end
    @test occursin(raw"\omega_{0}", string(latexify(w0)))
    @test !occursin(raw"\_{", string(latexify(w0)))
    @test occursin(raw"v_{x}\left( x, t \right)", string(latexify(vx)))
    @test !occursin(raw"\_{", string(latexify(vx)))
    for s in (string(latexify(eq)), repr(MIME"text/latex"(), eq))
        @test occursin(raw"v_{x}", s)
        @test occursin(raw"\omega_{0}", s)
        @test !occursin(raw"\_{", s)
    end
    # Powers of function variables must stay parenthesized
    @test occursin(raw"\left( v_{x}\left( x, t \right) \right)^{2}", string(latexify(vx^2)))
    # Unannotated multi-character names still go through the default Symbol path
    @variables plain_x
    @test occursin("\\mathtt{plain\\_x}", string(latexify(plain_x)))
end

# Custom-function arguments must keep outer recipe/caller Latexify options
@testset "latexwrapper argument formatting" begin
    @variables x y a[1:2] a_b
    @variables f(..) [latexwrapper = string]
    @test String(latexify(f(a[1]); env = :raw, index = :subscript)) ==
        raw"f\left( a_{1} \right)"
    @test String(latexify(f(1.23456789); env = :raw, fmt = "%.2f")) ==
        raw"f\left( 1.23 \right)"
    @test String(latexify(f(x * y); env = :raw, mult_symbol = raw"\times")) ==
        raw"f\left( x \times y \right)"
    @test String(latexify(f(a_b); env = :raw)) ==
        raw"f\left( \mathtt{a\_b} \right)"
    @test String(latexify(f(1.23456789e-9); env = :raw)) ==
        raw"f\left( 1.2346 \cdot 10^{-9} \right)"
end

# A derivative of a custom call multiplied by another factor must keep operand
# scope in the numerator, not as an unfenced differential operator times both.
@testset "latexwrapper derivative factor scope" begin
    @variables x
    @variables f(..) [latexwrapper = string]
    @variables g(..) [latexwrapper = string]
    D = Differential(x)
    @test String(latexify(D(f(x)) * g(x); env = :raw)) ==
        raw"\frac{\mathrm{d}f\left( x \right)}{\mathrm{d}x} ~ g\left( x \right)"
    @test String(latexify(f(x) * D(g(x)); env = :raw)) ==
        raw"\frac{\mathrm{d}g\left( x \right)}{\mathrm{d}x} ~ f\left( x \right)"
end
