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

@testset "substitute complex values (issue 1109)" begin
    @variables x y
    p = 0.4 + 1.7im * x

    r = substitute(p, Dict(x => 0.2 + 1.0im))
    @test r isa Complex{Num}
    @test unwrap_const(Symbolics.value(r)) ≈ -1.3 + 0.34im
    @test unwrap_const(Symbolics.value(real(r))) ≈ -1.3
    @test unwrap_const(Symbolics.value(imag(r))) ≈ 0.34

    @test substitute(p, Dict(x => 0.2 + 1.0im); fold = Val(true)) ==
        Complex{Num}(Num(-1.2999999999999998), Num(0.34))

    # substituting a non-constant complex value must combine the real and
    # imaginary parts of both substituted parts
    r2 = substitute(p, Dict(x => im * y))
    @test r2 isa Complex{Num}
    @test isequal(unwrap(real(r2)), unwrap(0.4 - 1.7y))
    @test isequal(unwrap(imag(r2)), unwrap(Num(0)))

    r3 = substitute(p, Dict(x => y + 2.0im * a))
    @test r3 isa Complex{Num}
    @test unwrap_const(Symbolics.value(substitute(r3, Dict(y => 1.0, a => 1.0)))) ≈
        0.4 + 1.7im * (1.0 + 2.0im)

    # substituting a complex value into a real Num expression (issue 1813)
    z = substitute(x + 1, Dict(x => 1.0 + 2.0im))
    @test z isa Complex{Num}
    @test unwrap_const(Symbolics.value(z)) == 2.0 + 2.0im
    @test unwrap_const(Symbolics.value(real(z))) == 2.0
    @test unwrap_const(Symbolics.value(imag(z))) == 2.0

    @test substitute(x, Dict(x => 1.0 + 2.0im)) isa Complex{Num}
    @test substitute(x * im, Dict(x => im)) isa Complex{Num}
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
