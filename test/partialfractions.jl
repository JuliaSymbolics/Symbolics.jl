using Symbolics
import Symbolics: partial_frac_decomposition

@variables x

# https://en.neurochispas.com/algebra/4-types-of-partial-fractions-decomposition-with-examples/
@test_broken isequal(partial_frac_decomposition((3x-1) / (x^2 + x - 6), x), 2/(x+3) + 1/(x-2))
@test_broken isequal(partial_frac_decomposition((9x^2 + 34x + 14) / ((x+2)*(x^2 - x - 12)), x), 3/(x+2) - 1/(x+3) + 7/(x-4))

@test_broken isequal(partial_frac_decomposition((2x^2 + 29x - 1) / ((2x+1)*(x - 2)^2), x), -4/(2x+1) + 3/(x-2) + 11/(x-2)^2)
