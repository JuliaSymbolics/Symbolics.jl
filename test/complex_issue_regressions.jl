using Test
using Symbolics
using SymbolicUtils
using SpecialFunctions
using LinearAlgebra

const SN = Symbolics.SymbolicNumber

@testset "historical complex issue regressions" begin
    @testset "#118 #1674 #718 #883 #1830 elementary complex functions" begin
        @variables x::Real z::Complex
        @test z isa SN
        @test !(z isa Complex{Num})

        phase = exp(im * x)
        @test phase isa SN
        @test !(phase isa Complex{Num})
        @test SymbolicUtils.operation(Symbolics.unwrap(phase)) === exp

        for (f, expected_op) in ((exp, exp), (sin, sin), (cos, cos), (log, log), (sqrt, sqrt))
            ex = f(z)
            @test ex isa SN
            @test !(ex isa Complex{Num})
            @test SymbolicUtils.operation(Symbolics.unwrap(ex)) === expected_op
        end

        shifted = sqrt(-im + x)
        logged = log(-im + x)
        @test shifted isa SN
        @test logged isa SN
    end

    @testset "#232 derivative of exp(im*x)" begin
        @variables x::Real
        D = Differential(x)
        dex = expand_derivatives(D(exp(im * x)))
        f = build_function(dex, x; expression = Val(false))
        @test f(0.37) ≈ im * exp(0.37im)
    end

    @testset "#327 #921 variable discovery" begin
        @variables t::Real x::Real u::Real v::Real z::Complex
        ex = x + t * Complex(u, v) + z
        vars = Set(Symbolics.get_variables(ex))
        for q in (x, t, u, v, z)
            @test any(vv -> isequal(vv, Symbolics.unwrap(q)), vars)
        end
        @test length(vars) == 5
        @test Set(Symbolics.get_variables(exp(im * x))) == Set([Symbolics.unwrap(x)])
    end

    @testset "#534 #905 #1109 #1813 substitution and domains" begin
        @variables z::Number f::Real
        sub1 = substitute(im * z, Dict(z => im); fold = Val(true))
        @test Symbolics.value(sub1) == -1
        ex = 0.4 + 1.7im * z
        sub2 = substitute(ex, Dict(z => 0.2 + 1.0im); fold = Val(true))
        @test Symbolics.value(sub2) ≈ -1.3 + 0.34im
        transfer = z^2 + 2z + 1
        got = substitute(transfer, Dict(z => 2pi * f * im))
        @test got isa SN

        @variables x::Real
        widened = substitute(x + 1, Dict(x => 1 + 2im); fold = Val(true))
        @test Symbolics.value(widened) == 2 + 2im
        @test substitute(x + 1, Dict(x => 2.0); fold = Val(true)) == 3.0
    end

    @testset "#159 #311 #354 #1016 #1391 build_function and expand" begin
        @variables z::Complex x::Real
        fz = build_function(z^2 + im * z, z; expression = Val(false))
        @test fz(1 + 2im) ≈ (1 + 2im)^2 + im * (1 + 2im)

        arr = [z, exp(im * z), 1 + im * z]
        fout, fout! = build_function(arr, z; expression = Val(false))
        @test fout(0.2 + 0.7im) ≈ [0.2 + 0.7im, exp(im * (0.2 + 0.7im)), 1 + im * (0.2 + 0.7im)]

        @test expand((x + im)^2) isa SN

        D = Differential(x)
        de = expand_derivatives(D(1 + exp(im * x)))
        dfun = build_function(de, x; expression = Val(false))
        @test dfun(0.0) ≈ im

        @variables a::Real b::Real
        cart = complex(a, b)
        cfun = build_function(cart, a, b; expression = Val(false))
        @test cfun(2.0, -3.0) == 2.0 - 3.0im
    end

    @testset "#341 expansion and #1116 degree" begin
        @variables z::Complex x::Real
        @test expand((z + 1)^3) isa SN
        @test Symbolics.degree(im + x, x) == 1
    end

    @testset "#777 rational construction semantics" begin
        @variables z::Complex
        @test (z / z) isa Number
        @test (z + im) / (z + im) isa Number
        @test_throws Exception (z + im) // (z + im)
        @test (1 // 2) * z isa SN
    end

    @testset "#800 printing" begin
        @variables z::Complex x::Real
        s1 = sprint(show, z)
        s2 = sprint(show, exp(im * x))
        @test occursin("z", s1)
        @test !occursin("real(", s1)
        @test !occursin("imag(", s1)
        @test occursin("exp", s2)
    end

    @testset "#832 real/imag simplification" begin
        @variables r1::Real r2::Real i1::Real i2::Real
        x1 = r1 + i1 * im
        x2 = r2 + i2 * im
        got = simplify(real(x1 * x2); expand = true)
        expected = r1 * r2 - i1 * i2
        @test isequal(simplify(got - expected; expand = true), 0)
    end

    @testset "#861 complex symbolic LinearAlgebra" begin
        @variables omega0::Real Omega::Real Delta::Real
        M = [-omega0 2im * Omega 0; -2im * Omega -omega0 2im * Delta; 0 -2im * Delta -omega0]
        d = det(M)
        fd = build_function(d, omega0, Omega, Delta; expression = Val(false))
        vals = (1.2, 0.4, -0.7)
        Mn = [-vals[1] 2im * vals[2] 0; -2im * vals[2] -vals[1] 2im * vals[3]; 0 -2im * vals[3] -vals[1]]
        @test fd(vals...) ≈ det(Mn)
    end

    @testset "#884 compact atomic rational expression" begin
        @variables z::Complex
        ex = 1 / (1 - z^10)
        @test ex isa SN
        txt = sprint(show, ex)
        @test !occursin("real(z)", txt)
        @test !occursin("imag(z)", txt)
    end

    @testset "#894 sound complex differential expression" begin
        @variables x::Real z(x)::Complex
        D = Differential(x)
        ex = x * D(z) + z
        @test ex isa SN
        @test SymbolicUtils.symtype(Symbolics.unwrap(ex)) <: Number
        eq = D(z) ~ ex
        @test eq isa Equation
    end

    @testset "#1199 #1485 complex equations and linear solve" begin
        @variables x::Real
        eq = x + 3 + im ~ 0
        @test eq isa Equation
        f = build_function(eq.lhs, x; expression = Val(false))
        @test f(2.0) == 5 + im
        sol = solve_for(x + im, x)
        @test iszero(simplify(sol + im))
    end

    @testset "#1487 left division" begin
        @variables x::Complex y::Complex
        ex = x \ y
        f = build_function(ex, x, y; expression = Val(false))
        @test f(1 + 2im, 3 - im) ≈ ((1 + 2im) \ (3 - im))
    end

    @testset "#1661 heterogeneous numeric symbolic domains" begin
        @variables t::Real x::Number y::Complex z(t)::Real
        v = [t, x, y, z]
        @test eltype(v) == SN
        @test length(v) == 4
    end

    @testset "#1372 dot follows Julia Hermitian semantics" begin
        @variables mass::Real qsqu::Real
        q2 = [0, 0, -(mass^2 + qsqu) / sqrt(qsqu) / 2,
              -im * (mass^2 + qsqu) / sqrt(qsqu) / 2]
        got = simplify(dot(q2, q2))
        direct = simplify(sum(q2[i] * q2[i] for i in eachindex(q2)))
        gdot = build_function(got, mass, qsqu; expression = Val(false))
        gdirect = build_function(direct, mass, qsqu; expression = Val(false))
        m, q = 1.4, 2.3
        qnum = [0, 0, -(m^2 + q) / sqrt(q) / 2, -im * (m^2 + q) / sqrt(q) / 2]
        @test gdot(m, q) ≈ dot(qnum, qnum)
        @test gdirect(m, q) ≈ sum(v * v for v in qnum)
        @test !isapprox(gdot(m, q), gdirect(m, q))
    end

    @testset "#465 complex symbolic arrays scalarize" begin
        @variables (v::Complex)[1:2]
        got = scalarize([1 2; 3 4] * v)
        @test length(got) == 2
        @test isequal(got[1], v[1] + 2v[2])
        @test isequal(got[2], 3v[1] + 4v[2])
    end

    @testset "#645 symbolic arrays and complex substitution" begin
        @variables r[1:1, 1:1]::Real phi[1:1, 1:1]::Real
        H = hankelh1.(0, r)
        Q = im .* cos.(phi)
        y = scalarize((Q .* H)[1])
        sub = substitute(y, Dict(r[1, 1] => 2.0, phi[1, 1] => 0.3); fold = Val(true))
        @test Symbolics.value(sub) ≈ im * cos(0.3) * hankelh1(0, 2.0)
    end

    @testset "#577 do not implicitly split complex equations" begin
        @variables z::Complex w::Complex
        eq = z ~ w
        @test eq isa Equation
        @test !(eq isa AbstractArray)
        f = build_function(eq.lhs - eq.rhs, z, w; expression = Val(false))
        @test f(1 + 2im, 0.5 - im) ≈ 0.5 + 3im
    end

    @testset "#558 #1011 complex differentiation contract" begin
        @variables t::Real w(t)::Complex
        D = Differential(t)
        dw = expand_derivatives(D(w))
        @test !iszero(dw)
        @test isequal(expand_derivatives(D(conj(w))), conj(dw))
        @test isequal(expand_derivatives(D(real(w))), real(dw))
        @test isequal(expand_derivatives(D(imag(w))), imag(dw))

        @variables z::Complex
        Dz = Differential(z)
        @test expand_derivatives(Dz(z)) == 1
        # Non-holomorphic projections remain unevaluated rather than silently claiming zero.
        @test Symbolics.is_derivative(expand_derivatives(Dz(conj(z))))
        @test Symbolics.is_derivative(expand_derivatives(Dz(real(z))))
        @test Symbolics.is_derivative(expand_derivatives(Dz(imag(z))))
    end
end
