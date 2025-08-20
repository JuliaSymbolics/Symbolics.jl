using Test
using Symbolics
import Nemo, Groebner
import Symbolics: partial_frac_decomposition

@variables x

# https://en.neurochispas.com/algebra/4-types-of-partial-fractions-decomposition-with-examples/
@test isequal(partial_frac_decomposition((3x-1) / (x^2 + x - 6), x), 2/(x+3) + 1/(x-2))
@test isequal(partial_frac_decomposition((9x^2 + 34x + 14) / ((x+2)*(x^2 - x - 12)), x), expand(3/(x+2) + 7/(x-4) - 1/(x+3)))

# https://tutorial.math.lamar.edu/Problems/Alg/PartialFractions.aspx
# can't handle leading coefficients being not 1 in denominator
@test isequal(partial_frac_decomposition((17x-53)/(x^2 - 2x - 15), x), 4/(x-5) + 13/(x+3))
@test_broken isequal(partial_frac_decomposition((34-12x)/(3x^2 - 10x - 8), x), -9(3x+2) - 1/(x-4))
@test isequal(partial_frac_decomposition((125 + 4x - 9x^2)/((x-1)*(x+3)*(x+4)), x), expand(6/(x-1) - 8/(x+3) - 7/(x+4)))
@test isequal(partial_frac_decomposition((10x+35)/((x+4)^2), x), 10/(x+4) + -5/(x+4)^2)
@test_broken isequal(partial_frac_decomposition((6x+5)/((2x-1)^2), x), 3/(2x-1) + 8/(2x-1)^2)
@test isequal(partial_frac_decomposition((7x^2-17x+38)/((x+6)*(x-1)^2), x), 8/(x+6) + -1/(x-1) + 4/(x-1)^2)
@test_broken isequal(partial_frac_decomposition((4x^2 - 22x + 7)/((2x+3)*(x-2)^2), x), 4/(2x+3) + -3/(x-2)^2)
@test_broken isequal(partial_frac_decomposition((3x^2 + 7x + 28)/(x*(x^2 + x + 7)), x), 4/x + (3-x)/(x^2+x+7)) # irrational roots
@test isequal(partial_frac_decomposition((4x^3 + 16x + 7)/(x^2 + 4)^2, x), 4x/(x^2+4) + 7/(x^2+4)^2)

# check valid expressions
@test partial_frac_decomposition(sin(x), x) === nothing
@test partial_frac_decomposition(x^2/(x-1), x) === nothing
@test partial_frac_decomposition(1/(x^2 + 2), x) === nothing # irrational roots, should eventually be fixed