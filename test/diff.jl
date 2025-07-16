using Symbolics
using Test
using Symbolics: value

# Derivatives
@variables t σ ρ β
@variables x y z
@variables uu(t) uuˍt(t) v(t)[1:3]
D = Differential(t)
D2 = Differential(t)^2
Dx = Differential(x)

@test Symbol(D(D(uu))) === Symbol("uuˍtt(t)")
@test Symbol(D(uuˍt)) === Symbol(D(D(uu)))
@test Symbol(D(v[2])) === Symbol("getindex(vˍt(t), 2)")

test_equal(a, b) = @test isequal(simplify(a), simplify(b))

@testset "ZeroOperator handling" begin
    @test_throws ErrorException Differential(0.1)(x)
    @test_throws ErrorException Differential(1)(x)
    @test_throws ErrorException Differential(2)(2x)
end

#@test @macroexpand(@derivatives D'~t D2''~t) == @macroexpand(@derivatives (D'~t), (D2''~t))

@test isequal(expand_derivatives(D(t)), 1)
@test isequal(expand_derivatives(D(D(t))), 0)

dsin = D(sin(t))
@test isequal(expand_derivatives(dsin), cos(t))

dcsch = D(csch(t))
@test isequal(expand_derivatives(dcsch), simplify(-coth(t) * csch(t)))

@test isequal(expand_derivatives(D(-7)), 0)
@test isequal(expand_derivatives(D(sin(2t))), simplify(cos(2t) * 2))
@test isequal(expand_derivatives(D2(sin(t))), simplify(-sin(t)))
@test isequal(expand_derivatives(D2(sin(2t))), simplify(-sin(2t) * 4))
@test isequal(expand_derivatives(D2(t)), 0)
@test isequal(expand_derivatives(D2(5)), 0)

# Chain rule
dsinsin = D(sin(sin(t)))
test_equal(expand_derivatives(dsinsin), cos(sin(t))*cos(t))

d1 = D(sin(t)*t)
d2 = D(sin(t)*cos(t))
@test isequal(expand_derivatives(d1), simplify(t*cos(t)+sin(t)))
@test isequal(expand_derivatives(d2), cos(t)^2-sin(t)^2)

eqs = [σ*(y-x),
       x*(ρ-z)-y,
       x*y - β*z]
jac = Symbolics.jacobian(eqs, [x, y, z])
test_equal(jac[1,1], -1σ)
test_equal(jac[1,2], σ)
test_equal(jac[1,3], 0)
test_equal(jac[2,1],  ρ - z)
test_equal(jac[2,2], -1)
test_equal(jac[2,3], -1x)
test_equal(jac[3,1], y)
test_equal(jac[3,2], x)
test_equal(jac[3,3], -1β)

# issue #545
z = t + t^2
#test_equal(expand_derivatives(D(z)), 1 + t * 2)

z = t-2t
#test_equal(expand_derivatives(D(z)), -1)

# Variable dependence checking in differentiation
@variables a(t) b(a)
@test !isequal(D(b), 0)
@test isequal(expand_derivatives(D(t)), 1)
@test isequal(expand_derivatives(Dx(x)), 1)

@variables x(t) y(t) z(t)

@test isequal(expand_derivatives(D(x * y)), simplify(y*D(x) + x*D(y)))
@test isequal(expand_derivatives(D(x * y)), simplify(D(x)*y + x*D(y)))

@test isequal(expand_derivatives(D(2t)), 2)
@test isequal(expand_derivatives(D(2x)), 2D(x))
@test isequal(expand_derivatives(D(x^2)), simplify(2 * x * D(x)))

# n-ary * and +
# isequal(Symbolics.derivative(Term(*, [x, y, z*ρ]), 1), y*(z*ρ))
# isequal(Symbolics.derivative(Term(+, [x*y, y, z]), 1), 1)

@test iszero(expand_derivatives(D(42)))
@test all(iszero, Symbolics.gradient(42, [t, x, y, z]))
@test all(iszero, Symbolics.hessian(42, [t, x, y, z]))
@test isequal(Symbolics.jacobian([t, x, 42], [t, x]),
              Num[1  0
                  Differential(t)(x)           1
                  0  0])

# issue 252
@variables beta, alpha, delta
@variables x1, x2, x3

# expression
tmp = beta * (alpha * exp(x1) * x2 ^ (alpha - 1) + 1 - delta) / x3
# derivative w.r.t. x1 and x2
t1 = Symbolics.gradient(tmp, [x1, x2])
@test_nowarn Symbolics.gradient(t1[1], [beta])

@variables t k
@variables x(t)
D = Differential(k)
@test Symbolics.tosymbol(value(D(x))) === Symbol("xˍk(t)")

using Symbolics
@variables t x(t)
∂ₜ = Differential(t)
∂ₓ = Differential(x)
L = .5 * ∂ₜ(x)^2 - .5 * x^2
@test isequal(expand_derivatives(∂ₓ(L)), -1 * x)
@test isequal(expand_derivatives(Differential(x)(L) - ∂ₜ(Differential(∂ₜ(x))(L))), -1 * (∂ₜ(∂ₜ(x)) + x))
@test isequal(expand_derivatives(Differential(x)(L) - ∂ₜ(Differential(∂ₜ(x))(L))), (-1 * x) - ∂ₜ(∂ₜ(x)))

@variables x2(t)
@test isequal(expand_derivatives(Differential(x)(2 * x + x2 * x)), 2 + x2)

@variables x y
@variables u(..)
Dy = Differential(y)
Dx = Differential(x)
dxyu = Dx(Dy(u(x,y)))
@test isequal(expand_derivatives(dxyu), dxyu)
dxxu = Dx(Dx(u(x,y)))
@test isequal(expand_derivatives(dxxu), dxxu)

using Symbolics, LinearAlgebra, SparseArrays
using Test

canonequal(a, b) = isequal(simplify(a), simplify(b))

# Calculus
@variables t σ ρ β
@variables x y z
@test isequal(
    (Differential(z) * Differential(y) * Differential(x))(t),
    Differential(z)(Differential(y)(Differential(x)(t)))
)

@test canonequal(
                 Symbolics.derivative(sin(cos(x)), x),
                 -sin(x) * cos(cos(x))
                )

Symbolics.@register_symbolic no_der(x)
@test canonequal(
                 Symbolics.derivative([sin(cos(x)), hypot(x, no_der(x))], x),
                 [
                  -sin(x) * cos(cos(x)),
                  x/hypot(x, no_der(x)) +
                      no_der(x)*Differential(x)(no_der(x))/hypot(x, no_der(x))
                 ]
                )

Symbolics.@register_symbolic intfun(x)::Int
@test Symbolics.symtype(intfun(x).val) === Int

eqs = [σ*(y-x),
       x*(ρ-z)-y,
       x*y - β*z]

∂ = Symbolics.jacobian(eqs,[x,y,z])
for i in 1:3
    ∇ = Symbolics.gradient(eqs[i],[x,y,z])
    @test canonequal(∂[i,:],∇)
end

@test all(canonequal.(Symbolics.gradient(eqs[1],[x,y,z]),[σ * -1,σ,0]))
@test all(canonequal.(Symbolics.hessian(eqs[1],[x,y,z]),0))

du = [x^2, y^3, x^4, sin(y), x+y, x+z^2, z+x, x+y^2+sin(z)]
reference_jac = sparse(Symbolics.jacobian(du, [x,y,z]))

@test findnz(Symbolics.jacobian_sparsity(du, [x,y,z]))[[1,2]] == findnz(reference_jac)[[1,2]]

let
    @variables t x(t) y(t) z(t)
    @test Symbolics.exprs_occur_in([x,y,z], x^2*y) == [true, true, false]
end

@test isequal(Symbolics.sparsejacobian(du, [x,y,z]), reference_jac)


function f!(res,u)
    (x,y,z)=u
    res.=[x^2, y^3, x^4, sin(y), x+y, x+z^2, z+x, x+y^2+sin(z)]
end
function f1!(res,u,a,b;c)
    (x,y,z)=u
    res.=[a*x^2, y^3, b*x^4, sin(y), c*x+y, x+z^2, a*z+x, x+y^2+sin(z)]
end

input=rand(3)
output=rand(8)

findnz(Symbolics.jacobian_sparsity(f!, output, input))[[1,2]] == findnz(reference_jac)[[1,2]]
findnz(Symbolics.jacobian_sparsity(f1!, output, input,1,2,c=3))[[1,2]] == findnz(reference_jac)[[1,2]]

input = rand(2,2)
function f2!(res,u,a,b,c)
    (x,y,z)=u[1,1],u[2,1],u[3,1]
    res.=[a*x^2, y^3, b*x^4, sin(y), c*x+y, x+z^2, a*z+x, x+y^2+sin(z)]
end

findnz(Symbolics.jacobian_sparsity(f!, output, input))[[1,2]] == findnz(reference_jac)[[1,2]]

# Check for failures due to du[4] undefined
function f_undef(du,u)
    du[1] = u[1]
    du[2] = u[2]
    du[3] = u[3] + u[4]
end
u0 = rand(4)
du0 = similar(u0)
sparsity_pattern = Symbolics.jacobian_sparsity(f_undef,du0,u0)
udef_ref = sparse([1 0 0 0
                    0 1 0 0
                    0 0 1 1
                    0 0 0 0])
findnz(sparsity_pattern)[[1,2]] == findnz(udef_ref)[[1,2]]

using Symbolics

rosenbrock(X) = sum(1:length(X)-1) do i
    100 * (X[i+1] - X[i]^2)^2 + (1 - X[i])^2
end
rosenbrock2(X,a;b) = sum(1:length(X)-1) do i
    a * (X[i+1] - X[i]^2)^2 + (b - X[i])^2
end

@variables a,b
X = [a,b]
input = rand(2)

spoly(x) = simplify(x, expand=true)
rr = rosenbrock(X)

reference_hes = Symbolics.hessian(rr, X)
@test findnz(sparse(reference_hes))[1:2] == findnz(Symbolics.hessian_sparsity(rr, X))[1:2]
@test findnz(sparse(reference_hes))[1:2] == findnz(Symbolics.hessian_sparsity(rosenbrock, input))[1:2]
@test findnz(sparse(reference_hes))[1:2] == findnz(Symbolics.hessian_sparsity(rosenbrock2, input,100,b=1))[1:2]

sp_hess = Symbolics.sparsehessian(rr, X)
@test findnz(sparse(reference_hes))[1:2] == findnz(sp_hess)[1:2]
@test isequal(map(spoly, findnz(sparse(reference_hes))[3]), map(spoly, findnz(sp_hess)[3]))

#96
@variables t x(t)[1:4] ẋ(t)[1:4]
expression = sin(x[1] + x[2] + x[3] + x[4]) |> Differential(t) |> expand_derivatives
expression2 = substitute(expression, Dict(collect(Differential(t).(x) .=> ẋ)))
@test isequal(expression2, (ẋ[1] + ẋ[2] + ẋ[3] + ẋ[4])*cos(x[1] + x[2] + x[3] + x[4]))

@test isequal(
    Symbolics.derivative(ifelse(signbit(b), b^2, sqrt(b)), b),
    ifelse(signbit(b), 2b,(SymbolicUtils.unstable_pow(2Symbolics.unwrap(sqrt(b)), -1)))
)

# Chain rule
#
let
    @variables t x(..) y(..) z(..)
    Dt = Differential(t)

    @test isequal(expand_derivatives(Dt(x(y(t)))),
                  Dt(y(t))*Differential(y(t))(x(y(t))))

    @test isequal(expand_derivatives(Dt(x(y(t), z(t)))),
                  Dt(y(t))*Differential(y(t))(x(y(t), z(t))) + Dt(z(t))*Differential(z(t))(x(y(t), z(t))))
end

@variables x y
@register_symbolic foo(x, y, z::Array)
D = Differential(x)
@test_throws ErrorException expand_derivatives(D(foo(x, y, [1.2]) * x^2))

@variables t x(t) y(t)
D = Differential(t)

eqs = [D(x) ~ x, D(y) ~ y + x]

sub_eqs = substitute(eqs, Dict([D(x)=>D(x), x=>1]))
@test sub_eqs == [D(x) ~ 1, D(y) ~ 1 + y]

@variables x y
@test substitute([x + y; x - y], Dict(x=>1, y=>2)) == [3, -1]


# 530#discussion_r825125589
let
    using Symbolics
    @variables u[1:2] y[1:1] t
    u = collect(u)
    y = collect(y)
    @test isequal(Symbolics.jacobian([u;u[1]^2; y], u), Num[1 0
                                                            0 1
                                                            2u[1] 0
                                                            0 0])
end

# make sure derivative(x[1](t), y) does not fail
let
    @variables t a(t)
    vars = collect(@variables(x(t)[1:1])[1])
    ps = collect(@variables(ps[1:1])[1])
    @test Symbolics.derivative(ps[1], vars[1]) == 0
    @test Symbolics.derivative(ps[1], a) == 0
    @test Symbolics.derivative(x[1], a) == 0
end

# 580
let
    @variables x[1:3]
    y = [sin.(x); cos.(x)]
    dj = Symbolics.jacobian(y, x)
    @test !iszero(dj)
    @test isequal(dj, Symbolics.jacobian(Symbolics.scalarize.(y), x))
    sj = Symbolics.sparsejacobian(y, x)
    @test !iszero(sj)
    @test isequal(sj, Symbolics.jacobian(Symbolics.scalarize.(y), x))
end

# substituting iv of differentials
@variables t t2 x(t)
D = Differential(t)
ex = D(x)
ex2 = substitute(ex, [t=>t2])
@test isequal(operation(Symbolics.unwrap(ex2)).x, t2)
ex3 = substitute(D(x) * 2 + x / t, [t=>t2])
xt2 = substitute(x, [t => t2])
@test isequal(ex3, xt2 / t2 + 2Differential(t2)(xt2))

# 581
#
let
    @variables x(t)[1:3]
    @test iszero(Symbolics.derivative(x[1], x[2]))
end

#908
#
let
    using Symbolics
    @variables t
    @test isequal(expand_derivatives(Differential(t)(im*t)), im)
    @test isequal(expand_derivatives(Differential(t)(t^2 + im*t)), 2t + im)
end

# 1262
#
let
    @variables t  b(t)
	D = Differential(t)
	expr = b - ((D(b))^2) * D(D(b))
	expr2 = D(expr)
	@test isequal(expand_derivatives(expr), expr)
    @test isequal(expand_derivatives(expr2), D(b) - (D(b)^2)*D(D(D(b))) - 2D(b)*(D(D(b))^2))
end

# 1126
#
let
    @syms y f(y) g(y) h(y)
    D = Differential(y)

    expr_gen = (fun) -> D(D(((-D(D(fun))) / g(y))))

    expr = expr_gen(g(y))
    # just make sure that no errors are thrown in the following, the results are to complicated to compare
    expand_derivatives(expr)
    expr = expr_gen(h(y))
    expand_derivatives(expr)

    expr = expr_gen(f(y))
    expand_derivatives(expr)
end

# Check `is_derivative` function
let
    @variables t X(t) Y(t)
    @syms a b
    D = Differential(t)
    my_f(x, y) = x^3 + 2y

    # Single expressions.
    @test !Symbolics.is_derivative(Symbolics.unwrap(D))
    @test !Symbolics.is_derivative(Symbolics.unwrap(t))
    @test !Symbolics.is_derivative(Symbolics.unwrap(X))
    @test !Symbolics.is_derivative(Symbolics.unwrap(a))
    @test !Symbolics.is_derivative(Symbolics.unwrap(1))

    # Composite expressions.
    @test Symbolics.is_derivative(Symbolics.unwrap(D(X)))
    @test !Symbolics.is_derivative(Symbolics.unwrap(D(X) + 3))
    @test Symbolics.is_derivative(Symbolics.unwrap(D(X + 2a*Y)))
    @test !Symbolics.is_derivative(Symbolics.unwrap(D(X) + D(Y)))
    @test !Symbolics.is_derivative(Symbolics.unwrap(my_f(X, D(Y))))
end

# Zeroth derivative (#1163)
let
    @variables t
    Dt = Differential(t)^0
    @test isequal(Dt, identity)
    test_equal(Dt(t + 2t^2), t + 2t^2)
end

# Check `Function` inputs for derivative (#1085)
let
    @variables x
    @testset for f in [sqrt, sin, acos, exp, cis]
        @test isequal(
            Symbolics.derivative(f, x),
            Symbolics.derivative(f(x), x)
        )
    end
end

# Check `Function` inputs throw for non-Num second input (#1085)
let
    @testset for f in [sqrt, sin, acos, exp, cis]
        @test_throws TypeError Symbolics.derivative(f, rand())
        @test_throws TypeError Symbolics.derivative(f, Val(rand(Int)))
    end
end

# Derivative of a `BasicSymbolic` (#1085)
let
    x = Symbolics.Sym{Int}(:x)
    @testset for f in [sqrt, sin, acos, exp]
        @test isequal(
            Symbolics.derivative(f, x),
            Symbolics.derivative(
                f,
                Symbolics.BasicSymbolic(x)
            )
        )
    end
end

# Check ssqrt, scbrt, slog
let
    @variables x
    D = Differential(x)
    @test isequal(expand_derivatives(D(Symbolics.ssqrt(1 + x ^ 2))), simplify((2x) / (2Symbolics.ssqrt(1 + x^2))))
    @test isequal(expand_derivatives(D(Symbolics.scbrt(1 + x ^ 2))), simplify((2x) / (3Symbolics.scbrt(1 + x^2)^2)))
    @test isequal(expand_derivatives(D(Symbolics.slog(1 + x ^ 2))), simplify((2x) / (1 + x ^ 2)))
end

# Hessian sparsity involving unknown functions
let
    @variables x₁ x₂ p q[1:1]
    expr = 3x₁^2 + 4x₁ * x₂
    @test Matrix(Symbolics.hessian_sparsity(expr, [x₁, x₂])) == [true true; true false]

    expr = 3x₁^2 + 4x₁ * x₂ + p
    @test Matrix(Symbolics.hessian_sparsity(expr, [x₁, x₂])) == [true true; true false]

    # issue 643: example test2_num
    expr = 3x₁^2 + 4x₁ * x₂ + q[1]
    @test Matrix(Symbolics.hessian_sparsity(expr, [x₁, x₂])) == [true true; true false]

    # Custom function: By default assumed to be non-linear
    myexp(x) = exp(x)
    @register_symbolic myexp(x)
    expr = 3x₁^2 + 4x₁ * x₂ + myexp(p)
    @test Matrix(Symbolics.hessian_sparsity(expr, [x₁, x₂])) == [true true; true false]
    expr = 3x₁^2 + 4x₁ * x₂ + myexp(x₂)
    @test Matrix(Symbolics.hessian_sparsity(expr, [x₁, x₂])) == [true true; true true]

    mylogaddexp(x, y) = log(exp(x) + exp(y))
    @register_symbolic mylogaddexp(x, y)
    expr = 3x₁^2 + 4x₁ * x₂ + mylogaddexp(p, 2)
    @test Matrix(Symbolics.hessian_sparsity(expr, [x₁, x₂])) == [true true; true false]
    expr = 3x₁^2 + 4x₁ * x₂ + mylogaddexp(3, p)
    @test Matrix(Symbolics.hessian_sparsity(expr, [x₁, x₂])) == [true true; true false]
    expr = 3x₁^2 + 4x₁ * x₂ + mylogaddexp(p, 2)
    @test Matrix(Symbolics.hessian_sparsity(expr, [x₁, x₂])) == [true true; true false]
    expr = 3x₁^2 + 4x₁ * x₂ + mylogaddexp(p, q[1])
    @test Matrix(Symbolics.hessian_sparsity(expr, [x₁, x₂])) == [true true; true false]
    expr = 3x₁^2 + 4x₁ * x₂ + mylogaddexp(p, x₂)
    @test Matrix(Symbolics.hessian_sparsity(expr, [x₁, x₂])) == [true true; true true]
    expr = 3x₁^2 + 4x₁ * x₂ + mylogaddexp(x₂, 4)
    @test Matrix(Symbolics.hessian_sparsity(expr, [x₁, x₂])) == [true true; true true]

    # Custom linear function: Possible to extend `Symbolics.linearity_1`/`Symbolics.linearity_2`
    myidentity(x) = x
    @register_symbolic myidentity(x)
    Symbolics.linearity_1(::typeof(myidentity)) = true
    expr = 3x₁^2 + 4x₁ * x₂ + myidentity(p)
    @test Matrix(Symbolics.hessian_sparsity(expr, [x₁, x₂])) == [true true; true false]
    expr = 3x₁^2 + 4x₁ * x₂ + myidentity(q[1])
    @test Matrix(Symbolics.hessian_sparsity(expr, [x₁, x₂])) == [true true; true false]
    expr = 3x₁^2 + 4x₁ * x₂ + myidentity(x₂)
    @test Matrix(Symbolics.hessian_sparsity(expr, [x₁, x₂])) == [true true; true false]

    mymul1plog(x, y) = x * (1 + log(y))
    @register_symbolic mymul1plog(x, y)
    Symbolics.linearity_2(::typeof(mymul1plog)) = (true, false, false)
    expr = 3x₁^2 + 4x₁ * x₂ + mymul1plog(p, q[1])
    @test Matrix(Symbolics.hessian_sparsity(expr, [x₁, x₂])) == [true true; true false]
    expr = 3x₁^2 + 4x₁ * x₂ + mymul1plog(x₂, q[1])
    @test Matrix(Symbolics.hessian_sparsity(expr, [x₁, x₂])) == [true true; true false]
    expr = 3x₁^2 + 4x₁ * x₂ + mymul1plog(q[1], x₂)
    @test Matrix(Symbolics.hessian_sparsity(expr, [x₁, x₂])) == [true true; true true]
end

# issue #555
let
    # first example
    @variables p[1:1] x[1:1]
    p = collect(p)
    x = collect(x)
    @test collect(Symbolics.sparsehessian(p[1] * x[1], x)) == [0;;]
    @test isequal(collect(Symbolics.sparsehessian(p[1] * x[1]^2, x)), [2p[1];;])

    # second example
    @variables a[1:2]
    a = collect(a)
    ex = (a[1]+a[2])^2
    @test Symbolics.hessian(ex, [a[1]]) == [2;;]
    @test collect(Symbolics.sparsehessian(ex, [a[1]])) == [2;;]
    @test collect(Symbolics.sparsehessian(ex, a)) == fill(2, 2, 2)
end

# issue #847
let
    @variables x[1:2] y[1:2]
    x = Symbolics.scalarize(x)
    y = Symbolics.scalarize(y)

    z = (x[1] + x[2]) * (y[1] + y[2])
    @test Symbolics.islinear(z, x)
    @test Symbolics.isaffine(z, x)

    z = (x[1] + x[2])
    @test Symbolics.islinear(z, x)
    @test Symbolics.isaffine(z, x)
end

# issue #790
let
    c(x) = [sum(x) - 1]
    @variables xs[1:2] ys[1:1]
    w = Symbolics.scalarize(xs)
    v = Symbolics.scalarize(ys)
    expr = dot(v, c(w))
    @test !Symbolics.islinear(expr, w)
    @test Symbolics.isaffine(expr, w)
    @test collect(Symbolics.hessian_sparsity(expr, w)) == fill(false, 2, 2)
end

# issue #749
let
    @variables x y
    @register_symbolic Base.FastMath.exp_fast(x, y)
    expr = Base.FastMath.exp_fast(x, y)
    @test !Symbolics.islinear(expr, [x, y])
    @test !Symbolics.isaffine(expr, [x, y])
    @test collect(Symbolics.hessian_sparsity(expr, [x, y])) == fill(true, 2, 2)
end

let
    @variables x
    expr = ifelse(x <= 0, -x, x)
    @test Symbolics.hessian_sparsity(expr, [x]) == [false;;]
    @test iszero(Symbolics.simplify(only(Symbolics.hessian(expr, [x]))))

    expr = ifelse(x <= 0, -x^2, x^2)
    @test Symbolics.hessian_sparsity(expr, [x]) == [true;;]
    @test !iszero(Symbolics.simplify(only(Symbolics.hessian(expr, [x]))))

    expr = ifelse(x <= 0, -x^2, x)
    @test Symbolics.hessian_sparsity(expr, [x]) == [true;;]
    @test !iszero(Symbolics.simplify(only(Symbolics.hessian(expr, [x]))))
end

@testset "`throw_no_derivative`" begin
    @variables t x(t) y::Bool z(..)
    D = Differential(t)
    @register_symbolic f(x)

    for expr in [D(f(x)), D(f(x)) + x^2, ifelse(y, D(f(x)), x^2), D(z(f(x))), z(D(f(x)))]
        @test_throws Symbolics.DerivativeNotDefinedError expand_derivatives(expr; throw_no_derivative = true)
        @test_throws Symbolics.DerivativeNotDefinedError Symbolics.derivative(expr, t; throw_no_derivative = true)
    end

    arr = [f(x) + t, x * t]
    vars = [x, t]
    @test_throws Symbolics.DerivativeNotDefinedError Symbolics.derivative(arr, t; throw_no_derivative = true)
    @test_throws Symbolics.DerivativeNotDefinedError Symbolics.gradient(arr[1], vars; throw_no_derivative = true)
    @test_throws Symbolics.DerivativeNotDefinedError Symbolics.jacobian(arr, vars; throw_no_derivative = true)
    @test_throws Symbolics.DerivativeNotDefinedError Symbolics.sparsejacobian(arr, vars; throw_no_derivative = true)
    @test_throws Symbolics.DerivativeNotDefinedError Symbolics.hessian(arr[1], vars; throw_no_derivative = true)
    @test_throws Symbolics.DerivativeNotDefinedError Symbolics.sparsehessian(arr[1], vars; throw_no_derivative = true)

    @variables x(t)[1:2]
    @test isequal(expand_derivatives(D(x[1]); throw_no_derivative = true), D(x[1]))
end

# issue #1452
let
    @variables x y z
    f(x, y, z) = x * y^2 + z
    @test Symbolics.hessian_sparsity(f(x, y, z), [x, y, z]) == [false true false; true true false; false false false]
    g(x, y, z) = x * y^2 + z
    @register_symbolic g(x, y, z)
    @test Symbolics.hessian_sparsity(g(x, y, z), [x, y, z]) == fill(true, 3, 3)
end

# Chain rule for array variables (issue #1383)
let
    @variables t x(t) y1(x) y(x)[1:1] Y(x)[1:1,1:1]
    Dt, Dx = Differential.([t, x])
    @test isequal(expand_derivatives(Dt(y1)), Dt(x) * Dx(y1)) # scalar
    @test isequal(expand_derivatives(Dt(y[1])), Dt(x) * Dx(y[1])) # same for vector var
    @test isequal(expand_derivatives(Dt(Y[1,1])), Dt(x) * Dx(Y[1,1])) # same for matrix arr
end

@testset "Derivatives of Array Variables" begin
    @variables p[1:9]
    vp = collect(p)

    f = - 1.5log(5 + p[3]) - 1.5log(7 + p[1]) - 1.5log(8 - p[2]) - 1.5log(9 - p[4]) + 0.08p[4]*p[9] -1.5log(5 + p[3]) - 1.5log(7 + p[1]) - 1.5log(8 - p[2]) - 1.5log(9 - p[4]) + 0.08p[4]*p[9]

    @test iszero(Symbolics.unwrap.(Symbolics.gradient(f, vp) .- Symbolics.gradient(f, p)))
    @test iszero(Symbolics.unwrap.(Symbolics.hessian(f, vp) .- Symbolics.hessian(f, p)))
    @test iszero(Symbolics.unwrap.(Symbolics.jacobian([f], vp) .- Symbolics.jacobian([f], p)))
end