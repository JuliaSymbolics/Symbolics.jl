using Test
using Symbolics
using SymbolicUtils
using Latexify
using LinearAlgebra
using SparseArrays

const SN = Symbolics.SymbolicNumber

@testset "complex feature surfaces" begin
    @testset "#1391 cis construction, differentiation, and codegen" begin
        @variables x::Real

        c = cis(x)
        @test c isa SN

        fexpr = 1 + c
        f = build_function(fexpr, x; expression = Val(false))
        @test f(0.37) ≈ 1 + cis(0.37)

        D = Differential(x)
        dex = expand_derivatives(D(fexpr))
        df = build_function(dex, x; expression = Val(false))
        @test df(0.37) ≈ im * cis(0.37)
    end

    @testset "#1674 compact conjugated complex exponential differentiation" begin
        @variables y::Real
        D = Differential(y)
        phase = exp(im * y)
        @test phase isa SN
        dex = expand_derivatives(D(conj(phase)))
        @test !Symbolics.is_derivative(dex)
        f = build_function(dex, y; expression = Val(false))
        @test f(0.37) ≈ conj(im * exp(0.37im))
    end

    @testset "#389 #416 complex Latexify" begin
        eqn(s) = "\\begin{equation}\n" * s * "\n\\end{equation}\n"
        @syms sx::Real
        @test string(latexify(im * sx)) == eqn(raw"\mathit{i} ~ \mathtt{sx}")

        @variables x::Real z::Complex a::Real b::Real
        @test string(latexify(z)) == eqn("z")
        @test string(latexify(im * x)) == eqn(raw"\mathit{i} ~ x")
        @test string(latexify(exp(im * x))) == eqn(raw"e^{\mathit{i} ~ x}")
        @test string(latexify(Complex(a, b))) == eqn(raw"a + b ~ \mathit{i}")

        @test sprint(show, MIME"text/latex"(), z) == "\$\$ " * eqn("z") * " \$\$"
        @test sprint(show, MIME"text/latex"(), exp(im * x)) == "\$\$ " * eqn(raw"e^{\mathit{i} ~ x}") * " \$\$"
    end

    @testset "numeric codomains narrow and widen correctly" begin
        @variables z::Complex
        @test real(z) isa Num
        @test imag(z) isa Num
        @test abs(z) isa Num
        @test abs2(z) isa Num
        @test angle(z) isa Num
        @test conj(z) isa SN
    end

    @testset "explicit Cartesian representation remains interoperable" begin
        @variables a::Real b::Real z::Complex
        cart = Complex(a, b)
        value = 0.7 + 1.3im
        for fop in (exp, sin, cos, log, sqrt)
            ex = fop(cart)
            fn = build_function(ex, a, b; expression = Val(false))
            @test fn(real(value), imag(value)) ≈ fop(value)
        end

        mixed = z + cart
        @test mixed isa SN
        fmixed = build_function(mixed, z, a, b; expression = Val(false))
        @test fmixed(0.2 - 0.4im, real(value), imag(value)) ≈ 0.2 - 0.4im + value

        promoted = [z, cart]
        @test eltype(promoted) == SN
    end

    @testset "general numeric wrapper reaches high-level differentiation" begin
        @variables z::Complex w::Complex

        dz = Symbolics.derivative(z^2 + im * z, z)
        @test dz isa SN
        @test iszero(simplify(dz - (2z + im)))

        g = Symbolics.gradient(z * w + im * z, [z, w])
        @test length(g) == 2
        @test iszero(simplify(g[1] - (w + im)))
        @test iszero(simplify(g[2] - z))

        J = Symbolics.jacobian([z^2 + w, im * z + w^2], [z, w])
        @test size(J) == (2, 2)
        @test iszero(simplify(J[1, 1] - 2z))
        @test isone(simplify(J[1, 2]))
        @test iszero(simplify(J[2, 1] - im))
        @test iszero(simplify(J[2, 2] - 2w))

        Js = Symbolics.sparsejacobian([z^2 + w, im * z + w^2], [z, w])
        @test Js isa SparseMatrixCSC
        @test iszero(simplify(Js[1, 1] - 2z))
        @test iszero(simplify(Js[2, 1] - im))

        H = Symbolics.hessian(im * z^2 + z * w, [z, w])
        @test size(H) == (2, 2)
        @test iszero(simplify(H[1, 1] - 2im))
        @test isone(simplify(H[1, 2]))
        @test isone(simplify(H[2, 1]))
        @test iszero(simplify(H[2, 2]))

        Hs = Symbolics.sparsehessian(im * z^2 + z * w, [z, w])
        @test Hs isa SparseMatrixCSC
        @test iszero(simplify(Hs[1, 1] - 2im))
        @test isone(simplify(Hs[1, 2]))
        @test isone(simplify(Hs[2, 1]))
    end

    @testset "mixed-domain differentiation" begin
        @variables x::Real y::Real z::Complex
        same(a, b) = iszero(simplify(a - b; expand = true))

        @testset "complex expressions of real variables" begin
            J = Symbolics.jacobian([im * x], [x])
            Js = Symbolics.sparsejacobian([im * x], [x])
            @test J isa Matrix{SN}
            @test Js isa SparseMatrixCSC{SN}
            @test all(same.(J, [im;;]))
            @test all(same.(Js, J))

            H = Symbolics.hessian(im * x^2, [x])
            Hs = Symbolics.sparsehessian(im * x^2, [x])
            @test H isa Matrix{SN}
            @test Hs isa SparseMatrixCSC{SN}
            @test all(same.(H, [2im;;]))
            @test all(same.(Hs, H))

            @test same(Symbolics.derivative(im * x, x), im)
            @test all(same.(Symbolics.gradient(im * x * y, [x, y]), [im * y, im * x]))
        end

        @testset "real expressions of complex variables" begin
            @test all(iszero, Symbolics.jacobian([x^2], [z]))
            @test all(iszero, Symbolics.hessian(x^2, [z]))
        end

        @testset "mixed expressions and variables" begin
            J = Symbolics.jacobian([x * z, x^2 + im * z], [x, z])
            Js = Symbolics.sparsejacobian([x * z, x^2 + im * z], [x, z])
            @test all(same.(J, [z x; 2x im]))
            @test all(same.(Js, J))

            H = Symbolics.hessian(x^2 * z + im * x, [x, z])
            Hs = Symbolics.sparsehessian(x^2 * z + im * x, [x, z])
            Hl = Symbolics.sparsehessian(x^2 * z + im * x, [x, z]; full = false)
            @test all(same.(H, [2z 2x; 2x 0]))
            @test all(same.(Hs, H))
            @test all(same.(Hl, [2z 0; 2x 0]))
        end

        @testset "Cartesian expressions of real variables" begin
            c = Complex(x, y)^2
            @test all(same.(Symbolics.jacobian([c], [x, y]), [2x + 2im * y 2im * x - 2y]))
            @test all(same.(Symbolics.hessian(c, [x, y]), [2 2im; 2im -2]))
        end

        @testset "real inputs return Num containers" begin
            @test Symbolics.jacobian([x^2 * y], [x, y]) isa Matrix{Num}
            @test Symbolics.sparsejacobian([x^2 * y], [x, y]) isa SparseMatrixCSC{Num}
            @test Symbolics.hessian(x^2 * y, [x, y]) isa Matrix{Num}
            @test Symbolics.sparsehessian(x^2 * y, [x, y]) isa SparseMatrixCSC{Num}
            @test Symbolics.gradient(x^2 * y, [x, y]) isa Vector{Num}
            @test Symbolics.derivative(x^2 * y, x) isa Num
        end
    end

    @testset "general numeric wrapper reaches symbolic linear algebra" begin
        @variables z::Complex w::Complex
        M = [z 1; 1 w]

        @test lu(M; check = false) isa LinearAlgebra.LU
        @test iszero(
            simplify_fractions(
                expand(
                    det(M; laplace = false) - det(M; laplace = true)
                )
            )
        )

        Minv = inv(M; laplace = false)
        ident = simplify.(M * Minv)
        @test isone(ident[1, 1])
        @test iszero(ident[1, 2])
        @test iszero(ident[2, 1])
        @test isone(ident[2, 2])

        Mex = exp([z zero(z); zero(w) w])
        @test Mex isa Symbolics.Arr{SN, 2}
    end

    @testset "complex symbolic arrays preserve result domains" begin
        @variables (v::Complex)[1:2]
        nv = norm(v)
        @test nv isa Num
    end

    @testset "complex symbolic linear systems" begin
        @variables z::Complex w::Complex

        scalar_sol = symbolic_linear_solve(z + im ~ 0, z)
        @test iszero(simplify(scalar_sol + im))

        sols = symbolic_linear_solve([z + w ~ 1, z - w ~ im], [z, w])
        @test length(sols) == 2
        @test iszero(simplify(sols[1] - (1 + im) / 2))
        @test iszero(simplify(sols[2] - (1 - im) / 2))
    end

    @testset "semi-polynomial forms preserve complex coefficients" begin
        @variables x::Real y::Real
        expr = im * x + (1 + im) * y + 2

        A, c = semilinear_form([expr], [x, y])
        @test iszero(simplify(A[1, 1] - im))
        @test iszero(simplify(A[1, 2] - (1 + im)))
        @test iszero(simplify((A * [x, y] + c)[1] - expr))

        qexpr = im * x^2 + (1 + im) * x * y + 2y + 3
        Aq, Bq, v2, cq = semiquadratic_form([qexpr], [x, y])
        @test iszero(simplify((Aq * [x, y] + Bq * v2 + cq)[1] - qexpr))
    end
end
