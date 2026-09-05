using Symbolics
using SymbolicUtils
using Test
using Symbolics: value, unwrap, dstar_derivative, dstar_jacobian, jacobian

# Standard derivative
@variables x
D = Differential(x)

@test isequal(dstar_derivative(x, x), expand_derivatives(D(x)))
@test isequal(dstar_derivative(2x, x), expand_derivatives(D(2x)))
@test isequal(dstar_derivative(x^2, x), expand_derivatives(D(x^2)))
@test isequal(expand(dstar_derivative(sin(cos(x))*cos(cos(x)), x)), expand_derivatives(D(sin(cos(x))*cos(cos(x)))))
@test isequal(expand(dstar_derivative(cos(sin(exp(2x)) + cos(exp(2x))), x)), expand(expand_derivatives(D(cos(sin(exp(2x)) + cos(exp(2x)))))))
@test isequal(expand(dstar_derivative(2sin(x^4 - x) + 3cos(2x), x)), expand(expand_derivatives(D(2sin(x^4 - x) + 3cos(2x)))))
@test isequal(dstar_derivative(2x*exp(x) + 2*exp(x), x), expand(expand_derivatives(D(2x*exp(x) + 2*exp(x)))))

# Standard Jacobian
@variables x y z
Dx = Differential(x)
Dy = Differential(y)
Dz = Differential(z)

@test isequal(dstar_jacobian([x,y], [x,y]), jacobian([x,y], [x,y]))
@test isequal(dstar_jacobian([x,y,z], [x,y,z]), jacobian([x,y,z], [x,y,z]))
@test isequal(dstar_jacobian([x,y,z], [x]), jacobian([x,y,z], [x]))
@test isequal(dstar_jacobian([x], [x,y,z]), jacobian([x], [x,y,z]))
@test isequal(dstar_jacobian([x*y*z], [x,y,z]), jacobian([x*y*z], [x,y,z]))
@test isequal(dstar_jacobian([x*y*z, x+y+z, sqrt(x^2 + y^2 + z^2)], [x,y,z]), jacobian([x*y*z, x+y+z, sqrt(x^2 + y^2 + z^2)], [x,y,z]))
@test isequal(dstar_jacobian([(x^2+y^2)^2, (x^2+y^2)^2 * y], [x,y]), jacobian([(x^2+y^2)^2, (x^2+y^2)^2 * y], [x,y]))
@test isequal(expand.(dstar_jacobian([(x^2 + y^2)*y, (x^2+y^2)*x^2 + (x^2+y^2)*y^2], [x,y])), expand.(jacobian([(x^2 + y^2)*y, (x^2+y^2)*x^2 + (x^2+y^2)*y^2], [x,y])))

# Edge case Jacobian
@test isequal(dstar_jacobian([x], [x,x]), jacobian([x], [x,x]))
@test isequal(dstar_jacobian([x,x], [x]), jacobian([x,x], [x]))
@test isequal(dstar_jacobian([x,x], [x,x]), jacobian([x,x], [x,x]))

# Array symbolics
@variables z[1:3]
@test isequal(dstar_jacobian(z,z), jacobian(z,z))
@test isequal(dstar_jacobian(2z, z), jacobian(2z, z))

# based on use in FastDifferentiation.jl and the D* paper for testing
function spherical_harmonics(max_l::Integer, x, y, z)
    Pc = Dict{Tuple{Int,Int}, Any}()
    Cc = Dict{Int, Any}()
    Sc = Dict{Int, Any}()

    function P(l, m)
        get!(Pc, (l, m)) do
            if l == 0 && m == 0
                1.0
            elseif l == m
                (1 - 2m) * P(m - 1, m - 1)
            elseif l == m + 1
                (2m + 1) * z * P(m, m)
            else
                ((2l - 1) / (l - m)) * z * P(l - 1, m) - ((l + m - 1) / (l - m)) * P(l - 2, m)
            end
        end
    end

    function C(m)
        get!(Cc, m) do
            m == 0 ? 1 : x * S(m - 1) + y * C(m - 1)
        end
    end

    function S(m)
        get!(Sc, m) do
            m == 0 ? 0 : x * C(m - 1) - y * S(m - 1)
        end
    end

    factorial_approx(n) = sqrt(2π * n) * (n / ℯ * sqrt(n * sinh(1 / n) + 1 / (810 * n^6)))^n
    N(l, m) = m == 0 ? sqrt(2l + 1 / (4π)) : sqrt((2l + 1) / 2π * factorial_approx(l - m) / factorial_approx(l + m))
    Y(l, m) = m < 0 ? N(l, -m) * P(l, -m) * S(-m) : N(l, m) * P(l, m) * C(m)

    return Num[Num(Y(l, m)) for l in 0:max_l-1 for m in -l:l]
end

@variables sx sy sz
sh_vars = [sx, sy, sz]

# max_l=4 -> 16 expressions
sh4 = spherical_harmonics(4, sx, sy, sz)
@test isequal(dstar_jacobian(sh4, sh_vars), jacobian(sh4, sh_vars))

# max_l=5 -> 25 expressions
# (more dominator/postdominator sharing than max_l=4 exercises).
sh5 = spherical_harmonics(5, sx, sy, sz)
@test isequal(dstar_jacobian(sh5, sh_vars), jacobian(sh5, sh_vars))

sh13 = spherical_harmonics(13, sx, sy, sz)
# large enough that dstar_jacobian and jacobian accumulate floating point errors, so sub in vars and use isapprox
let
    dj = dstar_jacobian(sh13, sh_vars)
    j = jacobian(sh13, sh_vars)
    subs = Dict(sx => 0.3, sy => 0.5, sz => 0.7)
    dvals = Symbolics.value.(substitute.(dj, (subs,)))
    jvals = Symbolics.value.(substitute.(j, (subs,)))
    @test isapprox(Float64.(dvals), Float64.(jvals); rtol=1e-8)
end