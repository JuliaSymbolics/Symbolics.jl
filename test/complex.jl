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

@testset "Complex{Num} power" begin
    @variables x y

    w = Complex(x, y)
    @test w^2 isa Complex{Num}

    result = w^(1 / 3)
    @test result isa Complex{Num}

    result_rat = w^(1 // 3)
    @test result_rat isa Complex{Num}

    r = real(w^(1 / 3))
    i = imag(w^(1 / 3))
    r_fn = Symbolics.build_function(r, [x, y]; expression = Val{false})
    i_fn = Symbolics.build_function(i, [x, y]; expression = Val{false})

    @test r_fn([8.0, 0.0]) ≈ 2.0
    @test i_fn([8.0, 0.0]) ≈ 0.0 atol = 1.0e-10

    @test r_fn([0.0, 8.0]) ≈ sqrt(3) atol = 1.0e-10
    @test i_fn([0.0, 8.0]) ≈ 1.0 atol = 1.0e-10

    for p in (0.5, 1 // 2)
        powered = w^p
        rf = Symbolics.build_function(real(powered), [x, y]; expression = Val{false})
        imf = Symbolics.build_function(imag(powered), [x, y]; expression = Val{false})
        for v in (1.0e200, 1.0e-200)
            @test complex(rf([v, v]), imf([v, v])) ≈ complex(v, v)^p
        end
    end

    @variables p
    c = (1.0 + 2.0im)^p
    @test c isa Complex{Num}
    @test SymbolicUtils.shape(Symbolics.unwrap(c)) == UnitRange{Int}[]
    @test SymbolicUtils.shape(Symbolics.unwrap(real(c))) == UnitRange{Int}[]
    @test SymbolicUtils.shape(Symbolics.unwrap(imag(c))) == UnitRange{Int}[]
end

@testset "power and division inference" begin
    @test @inferred(a^2) isa Num
    @test @inferred(a^(1 // 2)) isa Num
    @test @inferred(z / z) isa Complex{Num}
    @test @inferred(inv(z)) isa Complex{Num}
end
