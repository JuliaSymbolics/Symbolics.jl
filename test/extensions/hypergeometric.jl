using Symbolics
using HypergeometricFunctions
using Test

@testset "HypergeometricFunctions Extension" begin
    @testset "Symbolic registration and evaluation" begin
        # Create symbolic variables with default values
        @variables a=1.0 b=2.0 c=3.0 z=0.5
        @variables a₁=1.0 a₂=2.0 a₃=3.0 b₁=4.0 b₂=5.0

        # Test _₁F₁: symbolic call should produce symbolic result
        expr1 = HypergeometricFunctions._₁F₁(a, b, z)
        @test Symbolics.iscall(Symbolics.unwrap(expr1))

        # Substituting values should give same result as direct numerical call
        val1 = Symbolics.value(Symbolics.substitute(expr1, Dict(a => 1.0, b => 2.0, z => 0.5); fold = Val(true)))
        ref1 = HypergeometricFunctions._₁F₁(1.0, 2.0, 0.5)
        @test val1 ≈ ref1

        # Test _₂F₁
        expr2 = HypergeometricFunctions._₂F₁(a, b, c, z)
        @test Symbolics.iscall(Symbolics.unwrap(expr2))

        val2 = Symbolics.value(Symbolics.substitute(expr2, Dict(a => 1.0, b => 2.0, c => 3.0, z => 0.5); fold = Val(true)))
        ref2 = HypergeometricFunctions._₂F₁(1.0, 2.0, 3.0, 0.5)
        @test val2 ≈ ref2

        # Test _₃F₂
        expr3 = HypergeometricFunctions._₃F₂(a₁, a₂, a₃, b₁, b₂, z)
        @test Symbolics.iscall(Symbolics.unwrap(expr3))

        val3 = Symbolics.value(Symbolics.substitute(expr3, Dict(a₁ => 1.0, a₂ => 2.0, a₃ => 3.0, b₁ => 4.0, b₂ => 5.0, z => 0.5); fold = Val(true)))
        ref3 = HypergeometricFunctions._₃F₂(1.0, 2.0, 3.0, 4.0, 5.0, 0.5)
        @test val3 ≈ ref3
    end

    @testset "Derivative rules" begin
        @variables a=1.0 b=2.0 c=3.0 z=0.5
        @variables a₁=1.0 a₂=2.0 a₃=3.0 b₁=4.0 b₂=5.0
        D = Differential(z)

        # Test d/dz _₁F₁(a, b, z) = (a/b) * _₁F₁(a+1, b+1, z)
        expr1 = HypergeometricFunctions._₁F₁(a, b, z)
        deriv1 = expand_derivatives(D(expr1))

        # Evaluate derivative symbolically and compare with expected formula
        deriv1_val = Symbolics.value(Symbolics.substitute(deriv1, Dict(a => 1.0, b => 2.0, z => 0.5); fold = Val(true)))
        expected1 = (1.0 / 2.0) * HypergeometricFunctions._₁F₁(2.0, 3.0, 0.5)
        @test deriv1_val ≈ expected1

        # Test d/dz _₂F₁(a, b, c, z) = (a*b/c) * _₂F₁(a+1, b+1, c+1, z)
        expr2 = HypergeometricFunctions._₂F₁(a, b, c, z)
        deriv2 = expand_derivatives(D(expr2))

        deriv2_val = Symbolics.value(Symbolics.substitute(deriv2, Dict(a => 1.0, b => 2.0, c => 3.0, z => 0.5); fold = Val(true)))
        expected2 = (1.0 * 2.0 / 3.0) * HypergeometricFunctions._₂F₁(2.0, 3.0, 4.0, 0.5)
        @test deriv2_val ≈ expected2

        # Test d/dz _₃F₂
        expr3 = HypergeometricFunctions._₃F₂(a₁, a₂, a₃, b₁, b₂, z)
        deriv3 = expand_derivatives(D(expr3))

        deriv3_val = Symbolics.value(Symbolics.substitute(deriv3, Dict(a₁ => 1.0, a₂ => 2.0, a₃ => 3.0, b₁ => 4.0, b₂ => 5.0, z => 0.5); fold = Val(true)))
        expected3 = (1.0 * 2.0 * 3.0) / (4.0 * 5.0) * HypergeometricFunctions._₃F₂(2.0, 3.0, 4.0, 5.0, 6.0, 0.5)
        @test deriv3_val ≈ expected3
    end

    @testset "Chain rule" begin
        @variables a=1.0 b=2.0 t=0.5
        Dt = Differential(t)

        # d/dt _₁F₁(a, b, t²) = (a/b) * _₁F₁(a+1, b+1, t²) * 2t
        expr = HypergeometricFunctions._₁F₁(a, b, t^2)
        deriv = expand_derivatives(Dt(expr))

        deriv_val = Symbolics.value(Symbolics.substitute(deriv, Dict(a => 1.0, b => 2.0, t => 0.5); fold = Val(true)))
        # At t=0.5: inner derivative * outer derivative = (a/b) * _₁F₁(a+1, b+1, 0.25) * 2*0.5
        expected = (1.0 / 2.0) * HypergeometricFunctions._₁F₁(2.0, 3.0, 0.25) * 1.0
        @test deriv_val ≈ expected
    end

    @testset "Higher-order derivatives" begin
        @variables a=1.0 b=2.0 z=0.5
        D = Differential(z)
        D2 = D^2

        expr = HypergeometricFunctions._₁F₁(a, b, z)
        deriv2 = expand_derivatives(D2(expr))

        # d²/dz² _₁F₁(a,b,z) = (a*(a+1)) / (b*(b+1)) * _₁F₁(a+2, b+2, z)
        deriv2_val = Symbolics.value(Symbolics.substitute(deriv2, Dict(a => 1.0, b => 2.0, z => 0.5); fold = Val(true)))
        expected = (1.0 * 2.0) / (2.0 * 3.0) * HypergeometricFunctions._₁F₁(3.0, 4.0, 0.5)
        @test deriv2_val ≈ expected
    end

    @testset "Composition with other functions" begin
        @variables a=1.0 b=2.0 z=0.5
        D = Differential(z)

        # Product rule: d/dz [sin(z) * _₁F₁(a, b, z)]
        expr = sin(z) * HypergeometricFunctions._₁F₁(a, b, z)
        deriv = expand_derivatives(D(expr))

        deriv_val = Symbolics.value(Symbolics.substitute(deriv, Dict(a => 1.0, b => 2.0, z => 0.5); fold = Val(true)))
        # cos(z) * _₁F₁ + sin(z) * (a/b) * _₁F₁(a+1, b+1, z)
        expected = cos(0.5) * HypergeometricFunctions._₁F₁(1.0, 2.0, 0.5) +
                   sin(0.5) * (1.0 / 2.0) * HypergeometricFunctions._₁F₁(2.0, 3.0, 0.5)
        @test deriv_val ≈ expected
    end
end
