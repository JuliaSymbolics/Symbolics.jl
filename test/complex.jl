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

@testset "complex ~ returns a marked SplitComplexEquation" begin
    @variables t x y
    @variables ψ(..) w::Complex
    Dx = Differential(t)

    # numeric, symbolic, and dependent-variable sides all produce the marked type
    p_num = x + im * y ~ 1 + 2im
    p_sym = im ~ x + y * im
    p_dep = Dx(ψ(t, 1)) ~ im * ψ(t, 1)
    p_z = z ~ w
    for p in (p_num, p_sym, p_dep, p_z)
        @test p isa SplitComplexEquation
        @test p isa AbstractVector{Equation}
        @test iscomplexsplit(p)
        @test length(p) == 2
        @test size(p) == (2,)
        @test eltype(p) == Equation
        @test p[1] isa Equation && p[2] isa Equation
        @test p[1] == p.real_eq && p[2] == p.imag_eq
        @test p.original isa Equation
    end
    @test isequal(p_num, Equation[x ~ 1, y ~ 2])
    @test isequal(p_dep, Equation[Dx(ψ(t, 1)) ~ 0, 0 ~ ψ(t, 1)])
    @test isequal(p_z, Equation[real(z) ~ real(w), imag(z) ~ imag(w)])
    @test p_dep.original == Equation(Dx(ψ(t, 1)), im * ψ(t, 1))
    @test p_num.original == Equation(x + im * y, 1 + 2im)

    # behaves like the Vector{Equation} it replaced
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

    # recombining the parts recovers the original complex equation
    recombined = (Num(p_z[1].lhs) + im * Num(p_z[2].lhs)) ~
        (Num(p_z[1].rhs) + im * Num(p_z[2].rhs))
    @test iscomplexsplit(recombined)
    @test recombined.original == p_z.original

    # user-written groupings and non-split results are not marked
    @test !iscomplexsplit(Equation[p_dep[1], p_dep[2]])
    @test !iscomplexsplit([x ~ 1, y ~ 2])
    @test !iscomplexsplit(x ~ y)
    @test !iscomplexsplit(p_dep[1])
    @test (x ~ 1 + 2im) isa Equation
    @test (ψ(t) ~ 2im) isa Equation
end
