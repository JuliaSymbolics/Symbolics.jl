using Symbolics
using Test
using Latexify
using ReferenceTests

@variables x y z u(x)
Dx = Differential(x)
Dy = Differential(y)

# issue 260
@test Symbolics._toexpr(3 * x / y) == :((3x) / y)
@test Symbolics._toexpr(3 * x^y) == :(3x^y)

@test_reference "latexify_refs/inverse.txt" latexify(x^-1)

@test_reference "latexify_refs/frac1.txt" latexify((z + x * y^-1) / sin(z))
@test_reference "latexify_refs/frac2.txt" latexify((3x - 7y * z^23) * (z - z^2) / x)

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

@test_reference "latexify_refs/stable_mul_ordering1.txt" latexify(x * y)
@test_reference "latexify_refs/stable_mul_ordering2.txt" latexify(y * x)

@test_reference "latexify_refs/equation1.txt" latexify(x ~ y + z)
@test_reference "latexify_refs/equation2.txt" latexify(x ~ Dx(y + z))

@test_reference "latexify_refs/equation_vec1.txt" latexify([x ~ y + z
                                                            y ~ x - 3z])
@test_reference "latexify_refs/equation_vec2.txt" latexify([Dx(u) ~ z
                                                            Dx(y) ~ y * x])

@test_reference "latexify_refs/complex1.txt" latexify(x^2 - y^2 + 2im * x * y)
@test_reference "latexify_refs/complex2.txt" latexify(3im * x)
@test_reference "latexify_refs/complex3.txt" latexify(1 - x + (1 + 2x) * im;
                                                      imaginary_unit = "\\mathbb{i}")
