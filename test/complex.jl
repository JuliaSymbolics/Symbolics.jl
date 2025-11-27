using Symbolics, Test
using SymbolicUtils: metadata, symbolic_type, ScalarSymbolic, symtype, iscall, operation, arguments
using Symbolics: unwrap
using SymbolicIndexingInterface: getname, hasname

@variables a b::Real z::Complex (Z::Complex)[1:10]

@testset "types" begin
    @test a isa Num
    @test b isa Num
    @test eltype(Z) <: Complex{Real}

    for x in [z, Z[1], z+a, z*a, z^2, z/z] # z/z is sus
        @test x isa Complex{Num}
        @test real(x) isa Num
        @test imag(x) isa Num
        @test conj(x) isa Complex{Num}
    end

    # issue #314
    bi = a+a*im
    bs = substitute(bi, (Dict(a=>1.0))) # returns 1.0 + im
    typeof(bs) # Complex{Num}
    bv = Symbolics.value.(bs)
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

@testset "symbolic_type" begin
    # Test that ComplexTerm has correct symbolic_type classification
    # Regression test for bug where ComplexTerm returned NotSymbolic()
    # causing namespace issues in ModelingToolkit

    @variables x y p

    # Create a ComplexTerm
    complex_expr = p * x * (1 + im * y)
    ct = unwrap(complex_expr)

    # Verify it's a ComplexTerm
    @test ct isa Symbolics.ComplexTerm

    # Test symbolic_type returns ScalarSymbolic, not NotSymbolic
    @test symbolic_type(ct) == ScalarSymbolic()

    # Verify symtype is still correct
    @test Symbolics.symtype(ct) == Complex{Real}

    # Test that ComplexTerm is recognized as symbolic (needed for proper traversal)
    @test iscall(ct)
    @test operation(ct) == Complex{Real}
    @test length(arguments(ct)) == 2

    # Test with array parts (ComplexTerm can have array .re and .im)
    @variables x_arr[1:3] y_arr[1:3]
    ct_with_arrays = Symbolics.ComplexTerm{Real}(x_arr, y_arr)

    # ComplexTerm itself is still scalar even with array parts
    @test symbolic_type(ct_with_arrays) == ScalarSymbolic()
    @test Symbolics.symtype(ct_with_arrays) == Complex{Real}
end
