using Symbolics
using Test
using Latexify
using ReferenceTests

using DomainSets: Interval

@variables x y z u(x) dx h[1:10,1:10] hh(x,y)[1:10,1:10] gg(x,y)[1:10,1:10] [latexwrapper = string]
@variables AA(x) [latexwrapper = string] X₁(x) [latexwrapper = string]
@variables a[1:10]
Dx = Differential(x)
Dy = Differential(y)

# issue 260
@test Symbolics._toexpr(3*x/y) == :((3x) / y)
@test Symbolics._toexpr(3*x^y) == :(3x^y)

@test_reference "latexify_refs/inverse.txt" latexify(x^-1)

@test_reference "latexify_refs/integral1.txt" latexify(Integral(dx in Interval(0,1))(x))
@test_reference "latexify_refs/integral2.txt" latexify(Integral(dx in Interval(-Inf,Inf))(u^2))
@test_reference "latexify_refs/integral3.txt" latexify(Integral(dx in Interval(-z,u))(x^2))

@test_reference "latexify_refs/frac1.txt" latexify((z + x*y^-1) / sin(z))
@test_reference "latexify_refs/frac2.txt"  latexify((3x - 7y*z^23) * (z - z^2) / x)

@test_reference "latexify_refs/minus1.txt" latexify(x - y)
@test_reference "latexify_refs/minus2.txt" latexify(x - y * z)
@test_reference "latexify_refs/minus3.txt" latexify(sin(x+y-z))

@test_reference "latexify_refs/unary_minus1.txt" latexify(-y)
@test_reference "latexify_refs/unary_minus2.txt" latexify(-y * z)

@test_reference "latexify_refs/derivative1.txt" latexify(Dx(y))
@test_reference "latexify_refs/derivative2.txt" latexify(Dx(u))
@test_reference "latexify_refs/derivative3.txt" latexify(Dx(x^2 + y^2 + z^2))
@test_reference "latexify_refs/derivative4.txt" latexify(Dy(u))
@test_reference "latexify_refs/derivative5.txt" latexify(Dx(Dy(Dx(y))))

@test_reference "latexify_refs/stable_mul_ordering1.txt" latexify(x * y)
@test_reference "latexify_refs/stable_mul_ordering2.txt" latexify(y * x)

@test_reference "latexify_refs/equation1.txt" latexify(x ~ y + z)
@test_reference "latexify_refs/equation2.txt" latexify(x ~ Dx(y + z))
@test_reference "latexify_refs/equation5.txt" latexify(AA^2 + AA + 1 + X₁)

@test_reference "latexify_refs/equation_vec1.txt" latexify([
    x ~   y +  z
    y ~   x - 3z
])
@test_reference "latexify_refs/equation_vec2.txt" latexify([
    Dx(u) ~   z
    Dx(y) ~   y*x
])

@test_reference "latexify_refs/complex1.txt" latexify(x^2-y^2+2im*x*y)
@test_reference "latexify_refs/complex2.txt" latexify(3im*x)
@test_reference "latexify_refs/complex3.txt" latexify(1 - x + (1+2x)*im; imaginary_unit="\\mathbb{i}")
@test_reference "latexify_refs/complex4.txt" latexify(im * Symbolics.Term(sqrt, [2]))

@syms c
@test_reference "latexify_refs/complex5.txt" latexify((3+im/im)c)

@test_reference "latexify_refs/indices1.txt" latexify(h[10,10])
@test_reference "latexify_refs/indices2.txt" latexify(h[10,10], index=:bracket)

# test for https://github.com/JuliaSymbolics/Symbolics.jl/issues/1167
# note these tests need updating if/when https://github.com/korsbo/Latexify.jl/issues/331 is fixed
@test_reference "latexify_refs/indices3.txt" latexify(hh[10,10])
@test_reference "latexify_refs/indices4.txt" latexify(gg[10,10])

@test_reference "latexify_refs/indices5.txt" latexify(a'a)

@test !occursin("identity", latexify(Num(π))) # issue #1254
