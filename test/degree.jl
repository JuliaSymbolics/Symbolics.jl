using Symbolics
using Test

@variables x, y, z

@test isequal(degree(1), 0)
@test isequal(degree(x), 1)
@test isequal(degree(x, x), 1)
@test isequal(degree(x, y), 0)

@test isequal(degree(x*y, y), 1)
@test isequal(degree(x*y, x), 1)
@test isequal(degree(x*y, x*y), 1)
@test isequal(degree(x*y, z), 0)

@test isequal(degree(x*y^2+2*x+y^3*x^3), 6)
@test isequal(degree(x*y^2+2*x+y^3*x^3, y), 3)
@test isequal(degree(x*y^2+2*x+y^3*x^6, x), 6)

@test isequal(degree(x*y^2+2*x+y^3*x^6, z), 0)
@test isequal(degree(x*y^2+2*x+y^3*x^6+z, z), 1)
@test isequal(degree(x*y^2+2*x+y^3*x^6*z, z), 1)

@test isequal(degree(x^2 + y^2 +z^2, sin(z)), 0)
@test isequal(degree(sin(z), z), 0)
@test isequal(degree(sin(z)), 1)

@test isequal(degree(x^2*sin(z), x), 2)
@test isequal(degree(x+exp(z), x), 1)
