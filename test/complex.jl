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

@testset "substitute into Complex{Num} powers" begin
    @variables t
    w = 2 - 54cos(-t) + im*(-54sin(-t))
    @test w isa Complex{Num}
    # `Complex{Num} ^ literal` used to throw `UndefVarError: Pow`
    p = w^0.5
    @test p isa Complex{Num}
    u = 1 + 0.79p^(1/3) + 1/(0.79p^(1/3))
    # substituting used to throw an ArgumentError shape mismatch because
    # rebuilt `complex(re, im)` terms got shape `Unknown(-1)`
    r = substitute(u, Dict(t => 7))
    @test r isa Complex{Num}
    # with folding the result is fully numeric
    rv = unwrap_const(unwrap(substitute(u, Dict(t => 7); fold = Val(true))))
    wv = 2 - 54cos(-7.0) + im*(-54sin(-7.0))
    pv = wv^0.5
    uv = 1 + 0.79pv^(1/3) + 1/(0.79pv^(1/3))
    @test rv isa ComplexF64
    @test rv ≈ uv
end
