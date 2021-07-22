using Symbolics, Test

@variables a b::Real z::Complex{Real} (Z::Complex{Real})[1:10]

@testset "types" begin
    @test a isa Num
    @test b isa Num
    @test eltype(Z) <: Complex{Real}

    for x in [z, Z[1], z+a, z*a, z^2, z/z] # z/z is sus
        @test x isa Complex{Num}
        @test real(x) isa Num
        @test imag(x) isa Num
        @test conj(x) isa Complex{Num}
        @test sin(x) isa Complex{Num}
        @test sqrt(x) isa Complex{Num}
    end
end

@testset "repr" begin
    @test repr(z) == "z"
    @test repr(a + b*im) == "a + b*im"
end


@testset "derivatives" begin
    @variables x
    #issue 232
    @test isequal(Symbolics.derivative(exp(im*x),x), -sin(x)+im*cos(x))
end
