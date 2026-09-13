using Symbolics, Test
using SymbolicUtils: BasicSymbolic, metadata, symtype, unwrap_const
using Symbolics: unwrap
using SymbolicIndexingInterface: getname, hasname

const SN = Symbolics.SymbolicNumber

@variables a b::Real z::Complex (Z::Complex)[1:10]

@testset "atomic complex scalar representation" begin
    @test a isa Num
    @test b isa Num
    @test z isa SN
    @test Z[1] isa SN

    for x in (z, Z[1], z + a, z * a, z^2)
        @test x isa SN
        @test symtype(unwrap(x)) <: Number
        @test real(x) isa Num
        @test imag(x) isa Num
        @test conj(x) isa SN
    end
    @test isone(simplify(z / z))

    @test repr(z) == "z"
    phase = exp(im * a)
    @test phase isa SN
    @test SymbolicUtils.operation(unwrap(phase)) === exp
end

@testset "literal imaginary unit" begin
    ai = a * im
    ia = im * a

    @test ai isa SN
    @test ia isa SN
    @test !isdefined(Symbolics, :IM)
    @test Set(Symbolics.get_variables(ai)) == Set([unwrap(a)])
    @test Set(Symbolics.get_variables(ia)) == Set([unwrap(a)])

    @test Symbolics.value(substitute(ai, Dict(a => 2.0); fold = Val(true))) == 2.0im
    @test Symbolics.value(substitute(ia, Dict(a => 2.0); fold = Val(true))) == 2.0im

    f = build_function(ai, a; expression = Val(false))
    @test f(2.0) == 2.0im
end

@testset "elementary complex functions remain atomic" begin
    for f in (exp, sin, cos, log, sqrt)
        y = f(z)
        @test y isa SN
        @test !(y isa Complex{Num})
        @test SymbolicUtils.operation(unwrap(y)) === f
    end

    for f in (sqrt, log, exp)
        y = f(-im + a)
        @test y isa SN
        @test !(y isa Complex{Num})
        @test SymbolicUtils.operation(unwrap(y)) === f
    end
end

@testset "explicit Cartesian representation stays opt-in" begin
    cart = Complex(a, b)
    @test cart isa Complex{Num}
    @test isequal(real(cart), a)
    @test isequal(imag(cart), b)

    for ex in (
        im * a,
        a * im,
        a + 3im,
        3im + a,
        a - 3im,
        3im - a,
        a * (2 + 3im),
        (2 + 3im) * a,
        a / (2 + 3im),
        (2 + 3im) / a,
    )
        @test ex isa SN
        @test !(ex isa Complex{Num})
        @test symtype(unwrap(ex)) <: Number
    end
end

@testset "atomic substitution and metadata" begin
    bi = a + a * im
    @test bi isa SN

    bs = substitute(bi, Dict(a => 1.0); fold = Val(true))
    @test Symbolics.value(bs) == 1.0 + 1.0im

    @variables x::Complex
    @test !isnothing(metadata(unwrap(x)))
    @test_nowarn substitute(x + 1.0, x => 1.0im)
end

@testset "atomic naming" begin
    @variables t x::Complex y(t)::Complex q(a, b)::Complex
    @test hasname(x) && getname(x) == :x
    @test hasname(y) && getname(y) == :y
    @test hasname(q) && getname(q) == :q
    @test !hasname(2x)
    @test !hasname(x + y)
end

@testset "numeric wrapper lattice" begin
    r = SymbolicUtils.Sym{Symbolics.VartypeT}(:wrapper_real; type = Real)
    n = SymbolicUtils.Sym{Symbolics.VartypeT}(:wrapper_number; type = Number)
    c = SymbolicUtils.Sym{Symbolics.VartypeT}(:wrapper_complex; type = Complex{Real})

    @test Symbolics.wrapper_type(Real) === Num
    @test Symbolics.wrapper_type(Number) === SN
    @test Symbolics.wrapper_type(Complex{Real}) === SN

    @test Symbolics.wrap(r) isa Num
    @test Symbolics.wrap(n) isa SN
    @test Symbolics.wrap(c) isa SN

    wide_zero = zero(SN)
    @test wide_zero isa SN
    @test symtype(unwrap(wide_zero)) <: Real
    @test Symbolics.wrap(unwrap(wide_zero)) isa Num

    mixed = [Symbolics.wrap(r), Symbolics.wrap(c)]
    @test eltype(mixed) === SN
    @test mixed[1] isa SN
    @test symtype(unwrap(mixed[1])) <: Real
    @test symtype(unwrap(mixed[2])) <: Complex
end

abstract type ProbeNumericDomain <: Number end
@symbolic_wrap struct ProbeNumericWrapper <: ProbeNumericDomain
    val::BasicSymbolic{Symbolics.VartypeT}
end
SymbolicUtils.unwrap(x::ProbeNumericWrapper) = x.val

@testset "custom numeric wrapper specificity" begin
    p = SymbolicUtils.Sym{Symbolics.VartypeT}(:probe_numeric; type = ProbeNumericDomain)
    @test Symbolics.wrapper_type(ProbeNumericDomain) === ProbeNumericWrapper
    @test Symbolics.wrap(p) isa ProbeNumericWrapper
    @test Symbolics.unwrap(Symbolics.wrap(p)) === p
    @test Symbolics.wrapper_type(Real) === Num
    @test Symbolics.wrapper_type(Complex{Real}) === SN
end

@testset "atomic arrays and code generation" begin
    @variables z1::Complex z2::Complex
    A = [z1 z2; conj(z1) z1 + z2]
    @test eltype(A) <: Number

    f, f! = build_function(A, z1, z2; expression = Val(false))
    out = f(1.0 + 2.0im, 3.0 - 1.0im)
    expected = [1.0 + 2.0im 3.0 - 1.0im; 1.0 - 2.0im 4.0 + 1.0im]
    @test out == expected

    dest = similar(out)
    f!(dest, 1.0 + 2.0im, 3.0 - 1.0im)
    @test dest == expected
end
