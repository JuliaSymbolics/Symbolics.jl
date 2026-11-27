using Symbolics, Test, ForwardDiff
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

@testset "integer extremes" begin
    for x in (typemin(Int), typemax(Int), typemin(Int8), UInt(2), typemax(UInt), big(2)^80, 3 // 1, UInt(3) // UInt(1))
        @test rect(x) == tri(x) == 0
    end
    @test rect(UInt(0)) == tri(UInt(0)) == 1
    @test tri(UInt(1)) == 0
    @test rect(UInt(1) // UInt(2)) == 1 // 2
    @test heaviside(typemin(Int)) == 0
    @test heaviside(typemax(UInt)) == 1
    @test dirac(typemin(Int)) == dirac(typemax(UInt)) == 0
end

@testset "exact arguments near the kinks and their limits" begin
    m = typemax(Int)
    xs = (
        (m - 1) // m, -(m - 1) // m, m // (m - 1), 1 // m, -1 // m, (m ÷ 2) // m, m // 2,
        typemax(Int) // 1, (typemin(Int) + 1) // 1, typemin(Int) // 1, 1 // 0, -1 // 0,
        0 // 1, 1 // 1, -1 // 1, 1 // 2, -1 // 2, UInt(0) // UInt(1), typemax(UInt) // UInt(2),
        typemin(Int), typemax(Int), UInt(0), UInt(1), typemax(UInt), big(10)^100, true,
    )
    for x in xs, f in (heaviside, rect, tri)
        exact = big(x)
        @test f(x) ≈ f(exact isa Integer ? BigFloat(exact) : exact)
        root = rootfunction(f)(x)
        sides = (left_continuous_function(f)(x), right_continuous_function(f)(x))
        # the rebuilt piecewise form agrees away from the jumps
        (root == 0 && f !== tri) || @test sides[root < 0 ? 1 : 2] == f(x)
    end
    @test tri((m - 1) // m) === 1 // m
end

@testset "values at ±Inf" begin
    for x in (Inf, -Inf, Inf32, big(Inf))
        @test heaviside(x) == (x > 0)
        @test rect(x) == tri(x) == dirac(x) == Symbolics.dcosc(x) == 0
        @test isnan(diraccomb(x, 1))
        for f in (heaviside, rect, tri)
            side = rootfunction(f)(x) < 0 ? left_continuous_function(f) : right_continuous_function(f)
            @test side(x) == f(x)
        end
    end
end

@testset "second derivative of sinc at extreme arguments" begin
    dcosc = Symbolics.dcosc
    reference(x) = setprecision(BigFloat, 4096) do
        y = x isa Rational ? big(x.num) / big(x.den) : BigFloat(x)
        s, c = sinpi(y), cospi(y)
        (-big(π) * s - 2 * (c - s / big(π) / y) / y) / y
    end
    for x in (
            Int64(2)^53 + 1, -(Int64(2)^53 + 1), typemax(Int), typemin(Int), typemax(UInt),
            (Int64(2)^54 + 1) // 2, typemax(Int) // 2, (typemax(Int) - 1) // typemax(Int), 7 // 2,
            1.0e155, -1.0e155, floatmax(Float64), floatmax(Float32), floatmax(Float16),
        )
        expected = oftype(float(x), reference(x))
        @test dcosc(x) isa typeof(float(x))
        @test dcosc(x) ≈ expected rtol = 8 * eps(typeof(float(x)))
        @test signbit(dcosc(x)) == signbit(expected)
    end
    @test dcosc(Int64(2)^53 + 1) > 0
    @test dcosc(1.0e155) ≈ -2.0e-310
    @test isnan(dcosc(NaN))
end

@testset "diraccomb divisibility of exact arguments" begin
    m = typemax(Int)
    @test diraccomb(1 // m, 1 // (m - 1)) === 0.0
    @test_throws DomainError diraccomb(m // 1, 1 // m)
    @test diraccomb(Int64(2)^53 + 1, 2.0) === 0.0
    @test diraccomb(2.0^53, 3) === 0.0
    @test diraccomb(Float32(0.5), 1) === 0.0f0
    @test diraccomb(1 // 3, 1.0) === 0.0
    @test_throws DomainError diraccomb(2 // m, 1 // m)
    @test_throws DomainError diraccomb(Int64(2)^53 + 2, 2.0)
    @test_throws DomainError diraccomb(typemin(Int), UInt(2))
    @test_throws DomainError diraccomb(1.5, 1 // 2)
    @test isnan(diraccomb(1 // 0, 1))
    @test diraccomb(3, Inf) == 0
end

@testset "sinc's second derivative does not allocate for machine integers" begin
    allocated(x) = @allocated Symbolics.dcosc(x)
    # only types no wider than a word on any build: `Int128` is wider on 32-bit ones
    for x in (0, 3, -4, typemax(Int), typemin(Int), UInt(7), Int8(5))
        allocated(x)
        @test allocated(x) == 0
    end
end

@testset "derivative of tri near its kinks" begin
    f = build_function(d(tri(t)), t; expression = Val{false})
    @test f.([-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]) == [0, 0.5, 1, 0, -1, -0.5, 0]
    m = typemax(Int)
    for x in (2.0^-60, nextfloat(0.0), 1 // m, prevfloat(1.0), (m - 1) // m, UInt(1) // UInt(m))
        @test f(x) == -1
    end
    @test f(-2.0^-60) == f(-1 // m) == 1
    for x in (nextfloat(1.0), m // (m - 1), m, typemin(Int), typemax(UInt), UInt(2), Inf, -Inf)
        @test f(x) == 0
    end
    @test f(0.25f0) === -1.0f0
    @test isnan(f(NaN))
end

@testset "complex arguments of the sinc second derivative" begin
    reference(z) = setprecision(BigFloat, 512) do
        y = Complex{BigFloat}(z)
        s, c = sinpi(y), cospi(y)
        (-(big(π) * s) - 2 * (c - s / big(π) / y) / y) / y
    end
    f = build_function(d(d(sinc(t))), t; expression = Val{false})
    for z in (1.0 + 1.0im, 0.1 + 0.2im, -0.3im, 2.5 - 3.0im, 1 + 0im, 0.25f0 + 1.0f0im)
        @test f(z) isa typeof(float(z))
        @test f(z) ≈ reference(z) rtol = 16 * eps(typeof(real(float(z))))
    end
    @test f(0.0 + 0.0im) ≈ -π^2 / 3
    @test isnan(f(complex(NaN, 0.0)))
end

@testset "sinc's second derivative next to large integers and half-integers" begin
    reference(x) = setprecision(BigFloat, 4096) do
        y = big(x.num) / big(x.den)
        s, c = sinpi(y), cospi(y)
        (-(big(π) * s) - 2 * (c - s / big(π) / y) / y) / y
    end
    ε = 1 // big(10)^90
    for x in (big(10)^100 - ε, big(10)^100 + 1 + ε, big(10)^100 + 1 // 2 - ε, 1 - ε)
        @test Symbolics.dcosc(x) ≈ reference(x) rtol = 1.0e-60
    end
end

@testset "argument types are kept" begin
    for f in (heaviside, rect, tri, Symbolics.dcosc), T in (Float16, Float32, Float64, BigFloat)
        @test f(T(0.25)) isa T
        @test f(T(2)) isa T
        @test isnan(f(T(NaN)))
    end
    x = ForwardDiff.Dual(0.25, 1.0)
    @test ForwardDiff.partials(tri(x), 1) == -1
    @test ForwardDiff.value(rect(x)) == 1
end

@testset "second derivative of sinc" begin
    d2 = d(d(sinc(t)))
    @test isequal(d2, d(cosc(t)))
    f = build_function(d2, t; expression = Val{false})
    @test f(0.0) ≈ f(0) ≈ -π^2 / 3
    @test at(d2, t => 0) ≈ -π^2 / 3
    # series and closed form against a high-precision reference
    reference(x) = setprecision(BigFloat, 512) do
        y = big(x)
        -big(π)^2 * sinpi(y) / (big(π) * y) - 2 * (cospi(y) - sinpi(y) / (big(π) * y)) / y^2
    end
    for x in (1.0e-12, 1.0e-4, 0.1, -0.3, 0.4999, 0.5, 1.0, 7.25)
        @test f(x) ≈ reference(x) rtol = 1.0e-13
        @test Symbolics.dcosc(big(x)) ≈ reference(x) rtol = 1.0e-60
    end
    @test Symbolics.dcosc(Inf) == 0
    # no rule for the third derivative: unevaluated rather than wrong
    @test Symbolics.is_derivative(Symbolics.unwrap(d(d2)))
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
    @test isequal(d(tri(t)), Symbolics.dtri(t))
    @test isequal(d(d(tri(t))), dirac(t + 1) - 2dirac(t) + dirac(t - 1))
    for x in (-2.0, -1.5, -0.5, 0.25, 0.75, 1.5)
        @test at(d(tri(t)), t => x) ≈ (tri(x + 1.0e-6) - tri(x - 1.0e-6)) / 2.0e-6 atol = 1.0e-8
    end
    @test isequal(d(dirac(t)), dirac(t, 1))
    @test isequal(d(d(heaviside(t))), dirac(t, 1))
    @test isequal(d(dirac(t, 1)), dirac(t, 2))
    @test isequal(d(sinc(t)), cosc(t))
    @test isequal(d(sinc(a * t)), a * cosc(a * t))
    @test Symbolics.is_derivative(Symbolics.unwrap(d(diraccomb(t, T))))
    @test_throws OverflowError d(dirac(t, typemax(Int)))
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
    end
    @test rootfunction(heaviside)(0.0) == 0
    @test rootfunction(rect)(0.5) == rootfunction(rect)(-0.5) == 0
    @test all(iszero ∘ rootfunction(tri), (-1.0, 0.0, 1.0))
    for x in (Int64(2)^32, -Int64(2)^32, typemax(Int), typemin(Int), typemax(UInt), UInt(2))
        @test rootfunction(tri)(x) != 0
        @test sign(rootfunction(tri)(x)) == sign(x)
        @test rootfunction(rect)(x) > 0
        @test piecewise(tri, x) == tri(x) == 0
        @test piecewise(rect, x) == rect(x) == 0
    end
    # A root and two sides cannot express a point value: exactly at a jump the
    # reconstructed form takes a one-sided value instead of 1/2.
    @test piecewise(heaviside, 0.0) == 1
    @test heaviside(0.0) == 0.5
    @test piecewise(rect, 0.5) == piecewise(rect, -0.5) == 0
    @test rect(0.5) == rect(-0.5) == 0.5

    x = Symbolics.unwrap(t)
    for f in (heaviside, rect, tri)
        @test rootfunction(f)(x) isa Symbolics.SymbolicT
        @test left_continuous_function(f)(x) isa Symbolics.SymbolicT
        @test right_continuous_function(f)(x) isa Symbolics.SymbolicT
    end
end

# Both sides of an impulse are zero, so a registration would let a piecewise rewrite
# replace it by zero. They must stay opaque.
@testset "impulses are not registered discontinuities" begin
    for f in (dirac, diraccomb),
            g in (rootfunction, left_continuous_function, right_continuous_function)
        @test_throws MethodError g(f)
    end
end
