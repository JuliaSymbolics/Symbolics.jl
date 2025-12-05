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

# issue #1109: substituting complex values into expressions with complex coefficients
@testset "complex value substitution" begin
    @variables x y

    # Basic case from issue #1109
    p1 = 0.4 + 1.7im * x
    result1 = substitute(p1, Dict(x => 0.2 + 1.0im))
    expected1 = 0.4 + 1.7im * (0.2 + 1.0im)
    @test result1 isa Complex{Num}
    @test isapprox(unwrap_const(unwrap(real(result1))), real(expected1))
    @test isapprox(unwrap_const(unwrap(imag(result1))), imag(expected1))

    # Both real and imag parts have the variable
    p2 = x + 2.0im * x
    result2 = substitute(p2, Dict(x => 1.0 + 0.5im))
    expected2 = (1.0 + 0.5im) + 2.0im * (1.0 + 0.5im)
    @test isapprox(unwrap_const(unwrap(real(result2))), real(expected2))
    @test isapprox(unwrap_const(unwrap(imag(result2))), imag(expected2))

    # Real value substitution still works
    p3 = 0.4 + 1.7im * x
    result3 = substitute(p3, Dict(x => 0.5))
    expected3 = 0.4 + 1.7im * 0.5
    @test isapprox(unwrap_const(unwrap(real(result3))), real(expected3))
    @test isapprox(unwrap_const(unwrap(imag(result3))), imag(expected3))

    # Two variables with complex substitution
    p4 = x + y * im
    result4 = substitute(p4, Dict(x => 1.0 + 2.0im, y => 3.0 + 4.0im))
    expected4 = (1.0 + 2.0im) + (3.0 + 4.0im) * im
    @test isapprox(unwrap_const(unwrap(real(result4))), real(expected4))
    @test isapprox(unwrap_const(unwrap(imag(result4))), imag(expected4))

    # Symbolics.value should work on the result
    result_val = Symbolics.value(result1)
    @test result_val isa Complex
    @test isapprox(result_val, expected1)

    # simplify and expand should work
    @test_nowarn Symbolics.simplify(result1)
    @test_nowarn Symbolics.expand(result1)
end
