using Symbolics

# https://en.wikipedia.org/wiki/Taylor_series#List_of_Maclaurin_series_of_some_common_functions
@variables x
@test taylor(exp(x), x, 0:9) - sum(x^n//factorial(n) for n in 0:9) == 0
@test taylor(log(1-x), x, 0:9) - sum(-x^n/n for n in 1:9) == 0
@test taylor(log(1+x), x, 0:9) - sum((-1)^(n+1)*x^n/n for n in 1:9) == 0

@test taylor(1/(1-x), x, 0:9) - sum(x^n for n in 0:9) == 0
@test taylor(1/(1-x)^2, x, 0:8) - sum(n * x^(n-1) for n in 1:9) == 0
@test taylor(1/(1-x)^3, x, 0:7) - sum((n-1)*n*x^(n-2)/2 for n in 2:9) == 0
for α in (-1//2, 0, 1//2, 1, 2, 3)
    @test taylor((1+x)^α, x, 0:7) - sum(binomial(α, n)*x^n for n in 0:7) == 0
end

@test taylor(sin(x), x, 0:7) - sum((-1)^n/factorial(2*n+1) * x^(2*n+1) for n in 0:3) == 0
@test taylor(cos(x), x, 0:7) - sum((-1)^n/factorial(2*n) * x^(2*n) for n in 0:3) == 0
@test taylor(tan(x), x, 0:7) - taylor(taylor(sin(x), x, 0:7) / taylor(cos(x), x, 0:7), x, 0:7) == 0
@test taylor(asin(x), x, 0:7) - sum(factorial(2*n)/(4^n*factorial(n)^2*(2*n+1)) * x^(2*n+1) for n in 0:3) == 0
@test taylor(acos(x), x, 0:7) - (π/2 - taylor(asin(x), x, 0:7)) == 0 # TODO: make π/2 a proper fraction (like Num(π)/2)
@test taylor(atan(x), x, 0:7) - taylor(asin(x/√(1+x^2)), x, 0:7) == 0

@test taylor(sinh(x), x, 0:7) - sum(1/factorial(2*n+1) * x^(2*n+1) for n in 0:3) == 0
@test taylor(cosh(x), x, 0:7) - sum(1/factorial(2*n) * x^(2*n) for n in 0:3) == 0
@test taylor(tanh(x), x, 0:7) - (x - x^3/3 + 2/15*x^5 - 17/315*x^7) == 0

# around x ≠ 0
@test substitute(taylor(√(x), x, 1, 0:6), x => x + 1) - taylor(√(1+x), x, 0:6) == 0

# equations
eq = sin(2*x) ~ 2*sin(x)*cos(x)
eq = taylor(eq, x, 0:7)
eqs = taylor_coeff(eq, x) # should automatically expand to 7th order
@test length(eqs) == 7+1 && all(isequal(eq.lhs, eq.rhs) for eq in eqs)


# expand quintic equation around x=1
@variables ϵ
x_series = series(x, ϵ, 0:3)
x_coeffs = taylor_coeff(x_series, ϵ)
eq = x^5 + ϵ*x ~ 1
eqs = taylor_coeff(substitute(eq, x => x_series), ϵ, 0:3)
sol = x_coeffs .=> [1, -1//5, -1//25, -1//125] # e.g. https://ekamperi.github.io/mathematics/2020/06/21/perturbation-theory.html#MathJax-Element-39-Frame
eqs = substitute(eqs, Dict(sol))
@test all(isequal(eq.lhs, eq.rhs) for eq in eqs)
