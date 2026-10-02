using Test
using Symbolics

@testset "historical complex PR regressions" begin
    @testset "promotion preserves the real/wide wrapper lattice" begin
        @variables x::Real z::Number
        SN = Symbolics.SymbolicNumber

        @test promote_type(Float64, Num) === Num
        @test promote_type(ComplexF64, Num) === SN
        @test promote_type(Num, SN) === SN
        @test eltype([1.0, x]) === Num
        @test eltype([1.0 + 2.0im, x]) === SN
        @test z isa SN
    end

    @testset "#420 complex scalar cancellation" begin
        @variables x::Real
        @test isequal(simplify(x / im * im), x)
    end

    @testset "#908 #911 complex differentiation" begin
        @variables t::Real
        D = Differential(t)
        @test iszero(simplify(expand_derivatives(D(im * t)) - im))
        @test iszero(simplify(expand_derivatives(D(exp(im * t))) - im * exp(im * t)))
    end

    @testset "#1763 sinpi/cospi/sincospi" begin
        @variables x::Real z::Number
        SN = Symbolics.SymbolicNumber

        sx, cx = sincospi(x)
        @test sx isa Num
        @test cx isa Num
        @test sinpi(z) isa SN
        @test cospi(z) isa SN

        fs = build_function(sinpi(z), z; expression = Val(false))
        fc = build_function(cospi(z), z; expression = Val(false))
        v = 0.3 + 0.4im
        @test fs(v) ≈ sinpi(v)
        @test fc(v) ≈ cospi(v)
    end

    @testset "#1492 linear expansion stays generic" begin
        @variables x::Real y::Real

        a, b, islinear = Symbolics.linear_expansion(2im * x + im, x)
        @test islinear
        @test iszero(simplify(a - 2im))
        @test iszero(simplify(b - im))

        a, b, islinear = Symbolics.linear_expansion(im * x + im * y, x)
        @test islinear
        @test iszero(simplify(a - im))
        @test iszero(simplify(b - im * y))

        _, _, islinear = Symbolics.linear_expansion(im * x^2 + im * x, x)
        @test !islinear
    end

    @testset "#160 #919 #1326 complex build_function" begin
        @variables a::Real b::Real
        out = a + im * b
        f = build_function(out, (a, b); expression = Val(false))
        @test f((1.0, 2.0)) == 1.0 + 2.0im

        fruntime = build_function(1.0 + im * a, a; expression = Val(false))
        @test fruntime(1.0) == 1.0 + 1.0im

        fexpr = build_function(1 + im * a, a)
        fcompiled = eval(fexpr)
        @test Base.invokelatest(fcompiled, 1) == 1 + im
    end
end
