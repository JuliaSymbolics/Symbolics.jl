using Test
using Symbolics
using SymbolicUtils
using Nemo

const SN = Symbolics.SymbolicNumber

@testset "extended complex issue regressions" begin
    @testset "#622 #719 fractional powers do not fail during symbolic construction" begin
        @variables a::Real x::Real

        p1 = a^1.1
        p2 = (-a)^1.1
        @test p1 isa Num
        @test p2 isa Num

        # Historical #719 MWE: symbolic construction must not numerically probe the
        # negative coefficient and throw a DomainError while building the expression.
        p3 = (-csch(x)^2)^(1 // 2)
        @test p3 isa Number
    end

    @testset "#1741 atomic and explicit-Cartesian noninteger powers" begin
        @variables z::Complex p::Real x::Real y::Real

        for exponent in (1 / 3, 1 // 3)
            ex = z^exponent
            @test ex isa SN
            @test !(ex isa Complex{Num})
            f = build_function(ex, z; expression = Val(false))
            value = 0.7 + 1.3im
            @test f(value) ≈ value^exponent
        end

        # A concrete complex coefficient raised to a symbolic real power should remain
        # one atomic symbolic scalar rather than being decomposed into Cartesian slots.
        cp = (1.0 + 2.0im)^p
        @test cp isa SN
        cpf = build_function(cp, p; expression = Val(false))
        @test cpf(0.37) ≈ (1.0 + 2.0im)^0.37

        # Explicit Complex{Num} remains an opt-in Cartesian representation and must still
        # support the ordinary principal noninteger power semantics.
        w = Complex(x, y)
        @test w isa Complex{Num}
        for exponent in (1 / 3, 1 // 3)
            wp = w^exponent
            @test wp isa Complex{Num}
            wf = build_function(wp, x, y; expression = Val(false))
            value = 0.7 + 1.3im
            @test wf(real(value), imag(value)) ≈ value^exponent
        end

        # A symbolic real exponent preserves the explicit Cartesian representation and
        # uses the same principal polar branch.
        wsp = w^p
        @test wsp isa Complex{Num}
        wspf = build_function(wsp, x, y, p; expression = Val(false))
        value = 0.7 + 1.3im
        @test wspf(real(value), imag(value), 0.37) ≈ value^0.37
    end

    @testset "#1917 cubic/quartic solvers with non-rational coefficients" begin
        @variables x::Real

        cubic = x^3 + 0.5x + 1.0
        cubic_roots = eval.(Symbolics.toexpr.(symbolic_solve(cubic, x)))
        @test length(cubic_roots) == 3
        @test all(r -> abs(r^3 + 0.5r + 1.0) < 1e-10, cubic_roots)

        quartic = x^4 + 0.5x + 1.0
        quartic_roots = eval.(Symbolics.toexpr.(symbolic_solve(quartic, x)))
        @test length(quartic_roots) == 4
        @test all(r -> abs(r^4 + 0.5r + 1.0) < 1e-10, quartic_roots)

        @variables t::Real
        parametric = x^3 + x^2 * cos(t) + x * sin(t) + 1
        roots = Symbolics.get_roots_deg3(parametric, x)
        @test length(roots) == 3
        r = Symbolics.unwrap_const(Symbolics.value(substitute(Symbolics.wrap(roots[2]), Dict(t => 2); fold = Val(true))))
        @test abs(r^3 + r^2 * cos(2.0) + r * sin(2.0) + 1) < 1e-10
    end
end
