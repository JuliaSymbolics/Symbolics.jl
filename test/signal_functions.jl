using Symbolics, Test
using Symbolics: heaviside, dirac, rect, tri, diraccomb, rootfunction,
    left_continuous_function, right_continuous_function

@variables t a T q[1:3]
D = Differential(t)
d(ex) = expand_derivatives(D(ex))
at(ex, dict) = Symbolics.value(substitute(ex, dict; fold = Val(true)))

@testset "numeric values" begin
    @test heaviside.([-Inf, -1.0, -0.0, 0.0, 1.0, Inf]) == [0, 0, 0.5, 0.5, 1, 1]
    @test heaviside(0) == 0.5
    @test heaviside(1.0f0) isa Float32
    @test heaviside(0 // 1) === 1 // 2

    @test rect.([-1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0]) == [0, 0.5, 1, 1, 1, 0.5, 0]
    @test rect(prevfloat(0.5)) == rect(-prevfloat(0.5)) == 1
    @test rect(nextfloat(0.5)) == rect(-nextfloat(0.5)) == 0
    @test rect(1.0f0) isa Float32

    @test tri.([-2.0, -1.0, -0.5, 0.0, 0.25, 1.0, 2.0]) == [0, 0, 0.5, 1, 0.75, 0, 0]

    @test dirac(1.0) == dirac(-2) == dirac(1.0, 3) == 0
    @test_throws DomainError dirac(0)
    @test_throws DomainError dirac(-0.0)
    @test_throws DomainError dirac(0.0, 2)
    @test_throws DomainError dirac(1.0, -1)

    @test diraccomb(0.5, 1.0) == diraccomb(-2.5, 2) == 0
    @test_throws DomainError diraccomb(0.0, 1.0)
    @test_throws DomainError diraccomb(-6, 2)
    @test_throws DomainError diraccomb(1.0, -1.0)

    for f in (heaviside, rect, tri, dirac, Base.Fix2(diraccomb, 1.0))
        @test isnan(f(NaN))
    end
end

@testset "symbolic application" begin
    for f in (heaviside, rect, tri, dirac, sinc, cosc)
        ex = f(t)
        @test ex isa Num
        @test Symbolics.operation(Symbolics.unwrap(ex)) === f
        @test f(Symbolics.unwrap(t)) isa Symbolics.SymbolicT
        @test isequal(f.(collect(q)), [f(q[i]) for i in 1:3])
    end
    @test Symbolics.operation(Symbolics.unwrap(dirac(t, 2))) === dirac
    @test Symbolics.operation(Symbolics.unwrap(diraccomb(t, T))) === diraccomb
    @test diraccomb(t, 2) isa Num
    @test diraccomb(1.5, T) isa Num
    @test endswith(string(heaviside(t)), "heaviside(t)")
end

@testset "sinc(::Num)" begin
    # Previously traced into `sinpi(t) / (π * t)`, which evaluated to `NaN` at 0.
    @test string(sinc(t)) == "sinc(t)"
    @test at(sinc(t), t => 0) == 1
    @test at(sinc(t), t => 0.5) ≈ 2 / π
    @test at(d(sinc(t)), t => 0) == 0
    @test build_function(sinc(t), t; expression = Val{false})(0.0) == 1
    # complex arguments still go through Base
    @test sinc(Complex(t, t)) isa Complex{Num}
end

@testset "derivatives" begin
    @test isequal(d(heaviside(t)), dirac(t))
    @test isequal(d(heaviside(2t - a)), 2dirac(2t - a))
    @test isequal(d(rect(t)), dirac(t + 1 // 2) - dirac(t - 1 // 2))
    @test isequal(d(rect(a * t)), a * (dirac(a * t + 1 // 2) - dirac(a * t - 1 // 2)))
    @test isequal(d(tri(t)), rect(t + 1 // 2) - rect(t - 1 // 2))
    for x in (-2.0, -1.5, -0.5, 0.25, 0.75, 1.5)
        @test at(d(tri(t)), t => x) ≈ (tri(x + 1.0e-6) - tri(x - 1.0e-6)) / 2.0e-6 atol = 1.0e-8
    end
    @test isequal(d(dirac(t)), dirac(t, 1))
    @test isequal(d(d(heaviside(t))), dirac(t, 1))
    @test isequal(d(dirac(t, 1)), dirac(t, 2))
    @test isequal(d(sinc(t)), cosc(t))
    @test isequal(d(sinc(a * t)), a * cosc(a * t))
    @test at(d(d(sinc(t))), t => 0.3) ≈ (cosc(0.3 + 1.0e-6) - cosc(0.3 - 1.0e-6)) / 2.0e-6 rtol = 1.0e-6
    @test Symbolics.is_derivative(Symbolics.unwrap(d(diraccomb(t, T))))
end

@testset "substitute" begin
    @test at(heaviside(t), t => 0) == 0.5
    @test at(rect(t / a), Dict(t => 1, a => 2)) == 0.5
    @test at(tri(t), t => 0.25) == 0.75
    @test at(dirac(t - a), Dict(t => 1.0, a => 0.0)) == 0
    @test at(dirac(t, 1), t => 2.0) == 0
    @test at(diraccomb(t, T), Dict(t => 0.5, T => 1.0)) == 0
    @test_throws DomainError at(dirac(t), t => 0)
    # without folding, an impulse at 0 stays symbolic
    @test isequal(substitute(dirac(t), t => 0), dirac(Num(0)))
end

@testset "build_function" begin
    exs = [heaviside(t), rect(t), tri(t), sinc(t), cosc(t)]
    fns = (heaviside, rect, tri, sinc, cosc)
    for (ex, f) in zip(exs, fns)
        compiled = build_function(ex, t; expression = Val{false})
        evaled = eval(build_function(ex, t))
        for x in (-1.0, -0.5, 0.0, 0.3, 0.5, 2.0)
            @test compiled(x) == f(x)
            @test Base.invokelatest(evaled, x) == f(x)
        end
    end
    for ex in (dirac(t), dirac(t, 1), diraccomb(t, 1.0))
        compiled = build_function(ex, t; expression = Val{false})
        @test compiled(0.5) == 0
        @test_throws DomainError compiled(0.0)
    end
end

@testset "discontinuities" begin
    piecewise(f, args...) = rootfunction(f)(args...) < 0 ?
        left_continuous_function(f)(args...) : right_continuous_function(f)(args...)
    for x in -2.0:0.125:2.0
        x == 0 || @test piecewise(heaviside, x) == heaviside(x)
        abs(x) == 0.5 || @test piecewise(rect, x) == rect(x)
        @test piecewise(tri, x) == tri(x)
        x == 0 || @test piecewise(dirac, x) == dirac(x)
        x == 0 || @test piecewise(dirac, x, 1) == dirac(x, 1)
        isinteger(x) || @test piecewise(diraccomb, x, 1.0) == diraccomb(x, 1.0)
    end
    @test rootfunction(heaviside)(0.0) == 0
    @test rootfunction(rect)(0.5) == rootfunction(rect)(-0.5) == 0
    @test all(iszero ∘ rootfunction(tri), (-1.0, 0.0, 1.0))
    @test rootfunction(dirac)(0.0) == rootfunction(dirac)(0.0, 1) == 0
    @test all(iszero ∘ Base.Fix2(rootfunction(diraccomb), 2.0), (-2.0, 0.0, 2.0, 4.0))
    # a sign change at every impulse, so that events can locate it
    @test rootfunction(diraccomb)(1.9, 2.0) * rootfunction(diraccomb)(2.1, 2.0) < 0

    x = Symbolics.unwrap(t)
    for (f, args) in (
            (heaviside, (x,)), (rect, (x,)), (tri, (x,)), (dirac, (x,)),
            (dirac, (x, 1)), (diraccomb, (x, 1.0)),
        )
        @test rootfunction(f)(args...) isa Symbolics.SymbolicT
        @test left_continuous_function(f)(args...) isa Symbolics.SymbolicT
        @test right_continuous_function(f)(args...) isa Symbolics.SymbolicT
    end
end
