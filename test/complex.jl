using Symbolics, Test
using SymbolicUtils: metadata, unwrap_const
using Symbolics: unwrap
using SymbolicIndexingInterface: getname, hasname

@variables a b::Real z::Complex (Z::Complex)[1:10]

@testset "types" begin
    @test a isa Num
    @test b isa Num
    @test eltype(Z) <: Complex{Num}

    for x in [z, Z[1], z+a, z*a, z^2, z/z] # z/z is sus
        @test x isa Complex{Num}
        @test real(x) isa Num
        @test imag(x) isa Num
        @test conj(x) isa Complex{Num}
    end

    # issue #314
    bi = a+a*im
    bs = substitute(bi, (Dict(a=>1.0))) # returns 1.0 + im
    @test bs isa Complex{Num}
    bv = unwrap_const(Symbolics.value(bs))
    @test typeof(bv) == ComplexF64
end

@testset "repr" begin
    @test repr(z) == "z"
    @test repr(a + b*im) == "a + b*im"
end

@testset "metadata" begin
    z1 = z+1.0
    @test_nowarn substitute(z1, z=>1.0im)
    @test metadata(z1) == unwrap(z1.im).metadata
    @test metadata(z1) == unwrap(z1.re).metadata
    z2 = 1.0 + z*im
    @test isnothing(metadata(unwrap(z1.re)))
end

@testset "getname" begin
    @variables t a b x::Complex y(t)::Complex z(a, b)::Complex
    @test hasname(x)
    @test getname(x) == :x
    @test hasname(y)
    @test getname(y) == :y
    @test hasname(z)
    @test getname(z) == :z
    @test !hasname(2x)
    @test !hasname(x + y)
end

@testset "complex ~ still returns a plain Vector{Equation}" begin
    @variables x y
    p = x + im * y ~ 1 + 2im
    @test p isa Vector{Equation}
    @test isequal(p, Equation[x ~ 1, y ~ 2])
    push!(p, x ~ 3)
    @test length(p) == 3
    p[1] = y ~ 4
    @test isequal(p[1], y ~ 4)
end

@testset "split_complex_equation returns a SplitComplexEquation" begin
    @variables t x y
    @variables ψ(..) w::Complex
    Dx = Differential(t)

    p_num = split_complex_equation(x + im * y, 1 + 2im)
    p_sym = split_complex_equation(im, x + y * im)
    p_dep = split_complex_equation(Dx(ψ(t, 1)), im * ψ(t, 1))
    p_z = split_complex_equation(z, w)
    for p in (p_num, p_sym, p_dep, p_z)
        @test p isa SplitComplexEquation
        @test p isa AbstractVector{Equation}
        @test length(p) == 2
        @test size(p) == (2,)
        @test eltype(p) == Equation
        @test p[1] isa Equation && p[2] isa Equation
        @test p[1] == p.real_eq && p[2] == p.imag_eq
        @test p.original isa Equation
    end
    @test isequal(p_num, Equation[x ~ 1, y ~ 2])
    @test isequal(p_sym, Equation[0 ~ x, 1 ~ y])
    @test isequal(p_dep, Equation[Dx(ψ(t, 1)) ~ 0, 0 ~ ψ(t, 1)])
    @test isequal(p_z, Equation[real(z) ~ real(w), imag(z) ~ imag(w)])
    @test p_dep.original == Equation(Dx(ψ(t, 1)), im * ψ(t, 1))
    @test p_num.original == Equation(x + im * y, 1 + 2im)

    # the marked pair holds the same equations as ~
    for (a, b) in ((x + im * y, 1 + 2im), (im, x + y * im), (Dx(ψ(t, 1)), im * ψ(t, 1)), (z, w))
        @test isequal(collect(split_complex_equation(a, b)), a ~ b)
    end

    @test [e for e in p_dep] == [p_dep[1], p_dep[2]]
    @test collect(p_dep) isa Vector{Equation}
    @test p_dep[end] == p_dep[2]
    @test p_dep[1:2] isa Vector{Equation}
    @test p_dep[1:2] == [p_dep[1], p_dep[2]]
    @test_throws BoundsError p_dep[0]
    @test_throws BoundsError p_dep[3]
    @test first(p_dep) == p_dep[1]
    @test last(p_dep) == p_dep[2]
    @test isequal(map(eq -> eq.lhs, p_dep), [p_dep[1].lhs, p_dep[2].lhs])
    @test isequal(Symbolics.lhss(p_dep), [p_dep[1].lhs, p_dep[2].lhs])
    @test vcat(p_dep) isa Vector{Equation}
    @test vcat(p_dep) == [p_dep[1], p_dep[2]]
    @test vcat(p_dep, x ~ y) isa Vector{Equation}
    @test length(vcat(p_dep, x ~ y)) == 3
    @test reduce(vcat, [p_dep, p_sym]) isa Vector{Equation}
    @test reduce(vcat, [p_dep, p_sym]) == [p_dep[1], p_dep[2], p_sym[1], p_sym[2]]
    @test p_dep == collect(p_dep)
    @test isequal(p_dep, collect(p_dep))
    @test hash(p_dep) == hash(collect(p_dep))

    # user-written groupings and the output of ~ are not marked
    @test !(Equation[p_dep[1], p_dep[2]] isa SplitComplexEquation)
    @test !((Dx(ψ(t, 1)) ~ im * ψ(t, 1)) isa SplitComplexEquation)

    # without a split, the result is the same single Equation as ~
    @test isequal(split_complex_equation(x, 1 + 2im), x ~ 1 + 2im)
    @test split_complex_equation(x, 1 + 2im) isa Equation
    @test isequal(split_complex_equation(ψ(t), 2im), ψ(t) ~ 2im)
    @test isequal(split_complex_equation(x, y), x ~ y)
    @test_throws ErrorException split_complex_equation(1 + 2im, 3im)
end
