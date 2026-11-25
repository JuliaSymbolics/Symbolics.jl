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
    @variables x y n::Number
    p = x + im * y ~ 1 + 2im
    @test p isa Vector{Equation}
    @test isequal(p, Equation[x ~ 1, y ~ 2])
    push!(p, x ~ 3)
    @test length(p) == 3
    p[1] = y ~ 4
    @test isequal(p[1], y ~ 4)
    @test (complex(Num(n), y) ~ 0) == Equation[Num(n) ~ 0, y ~ 0]
    @test (0 ~ complex(Num(n), y)) == Equation[0 ~ Num(n), 0 ~ y]
end

@testset "split_complex_equation with Number-typed symbols" begin
    @variables n::Number x y
    original = Equation(n, x + im * y)
    wrapped = split_complex_equation(Num(n), x + im * y)
    unwrapped = split_complex_equation(original.lhs, original.rhs)
    from_equation = split_complex_equation(original)
    values = Dict{Any, Number}(unwrap(n) => 1 + 2im, unwrap(x) => 1, unwrap(y) => 2)

    for p in (wrapped, unwrapped, from_equation)
        @test p == Equation[real(original.lhs) ~ x, imag(original.lhs) ~ y]
        @test p.original == original
    end
    for eq in (original, wrapped..., unwrapped..., from_equation...)
        @test unwrap_const(substitute(eq.lhs, values; fold = Val(true))) ==
            unwrap_const(substitute(eq.rhs, values; fold = Val(true)))
    end
end

@testset "split_complex_equation with Number-typed complex components" begin
    @variables n::Number y
    z = complex(Num(n), y)
    for (lhs, rhs) in ((z, 0), (0, z))
        @testset "$lhs ~ $rhs" begin
            p = split_complex_equation(lhs, rhs)
            @test p == (lhs ~ rhs)
            @test split_complex_equation(p.original) == p
        end
    end
end

@testset "split_complex_equation requires scalar operands" begin
    @variables A[1:2] B[1:2] (Z::Complex)[1:2] (W::Complex)[1:3]
    for (lhs, rhs, eq) in ((A, B, A ~ B), (Z, W, Equation(Z, W)))
        @test_throws ArgumentError split_complex_equation(lhs, rhs)
        @test_throws ArgumentError split_complex_equation(eq.lhs, eq.rhs)
        @test_throws ArgumentError split_complex_equation(eq)
    end
    @test_throws ArgumentError split_complex_equation([1, 2], 0)
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

    # unwrapped sides and an Equation argument split like the wrapped sides
    for p in (p_num, p_sym, p_dep, p_z)
        @test isequal(split_complex_equation(p.original.lhs, p.original.rhs), p)
        @test isequal(split_complex_equation(p.original), p)
    end
    @test isequal(split_complex_equation(x ~ y), x ~ y)

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

    # copy keeps the marked pair; equality ignores .original
    @test copy(p_dep) isa SplitComplexEquation
    @test copy(p_dep) == p_dep
    @test p_dep == SplitComplexEquation(p_dep.real_eq, p_dep.imag_eq, x ~ y)
    @test hash(p_dep) == hash(SplitComplexEquation(p_dep.real_eq, p_dep.imag_eq, x ~ y))

    @test_throws MethodError push!(p_dep, x ~ 3)
    @test !((Dx(ψ(t, 1)) ~ im * ψ(t, 1)) isa SplitComplexEquation)

    @test split_complex_equation(x, 1 + 0im) == (x ~ 1)
    @test split_complex_equation(1 + 0im, x) == (1 ~ x)
    @test split_complex_equation(im * x, 2im) == (x ~ 2)
    @test split_complex_equation(2im, im * x) == (2 ~ x)
    @test isequal(split_complex_equation(x, y), x ~ y)
    @test_throws ErrorException split_complex_equation(1 + 2im, 3im)
end

@testset "split_complex_equation rejects unequal numeric parts" begin
    @variables x
    for (side, candidate) in ((x, 1), (im * x, 2))
        for (lhs, rhs) in ((side, 1 + 2im), (1 + 2im, side))
            values = Dict(x => candidate)
            @test unwrap_const(unwrap(substitute(lhs, values))) !=
                unwrap_const(unwrap(substitute(rhs, values)))
            @test_throws ArgumentError split_complex_equation(lhs, rhs)
        end
    end
end

@testset "split_complex_equation equal and unequal numeric parts" begin
    @variables t x
    @variables ψ(..)
    # equal NaN parts are accepted, matching ~
    @test split_complex_equation(NaN + im * x, NaN + 2im) == (x ~ 2.0)
    @test split_complex_equation(x + NaN * im, 2 + NaN * im) == (x ~ 2.0)
    @test split_complex_equation(-0.0 + im * x, 0.0 + 2im) == (x ~ 2.0)
    @test isequal(split_complex_equation(NaN + im * x, NaN + 2im), (NaN + im * x) ~ (NaN + 2im))
    # a NaN against a number stays unequal
    @test_throws ArgumentError split_complex_equation(NaN + im * x, 1 + 2im)
    @test_throws ArgumentError split_complex_equation(x + NaN * im, 2 + 1im)
    # Bool parts of a Complex{Bool} are reported as numbers
    err = try
        split_complex_equation(ψ(t, 0), im)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test err.msg == "Unequal numeric imaginary parts in equation: 0 and 1"
    err = try
        split_complex_equation(im, 1 + x * im)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test err.msg == "Unequal numeric real parts in equation: 0 and 1"
end

@testset "power and division inference" begin
    @test @inferred(a^2) isa Num
    @test @inferred(a^(1 // 2)) isa Num
    @test @inferred(z / z) isa Complex{Num}
    @test @inferred(inv(z)) isa Complex{Num}
end
