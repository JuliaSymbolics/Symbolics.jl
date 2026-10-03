using Symbolics: diff2term, value
using ModelingToolkit
using Test

# Test diff2term basic functionality (unit-specific tests removed as API changed in MTK v11)

@variables x y t
D = Differential(t)
# Test simple variable derivative
@test Symbolics.diff2term(D(x)) isa Num

@variables x(t)[1:2] p[1:2, 1:2]
@test ModelingToolkit.is_diff_equation(D(x) ~ p*x)
