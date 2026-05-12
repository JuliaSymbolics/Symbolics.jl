using Symbolics, NaNMath, Test

function discontinuity_eval(fn, args...)
    if rootfunction(fn)(args...) < 0
        left_continuous_function(fn)(args...)
    else
        right_continuous_function(fn)(args...)
    end
end

@testset "abs" begin
    for x in -1.0:0.001:1.0
        @test abs(x) ≈ discontinuity_eval(abs, x)
    end
end

@testset "$(nameof(f))" for f in (mod, rem, div)
    y = 0.7
    for x in -2y:0.001:2y
        @test f(x, y) ≈ discontinuity_eval(f, x, y)
    end
end

@testset "$(nameof(f))" for f in (min, max, NaNMath.min, NaNMath.max, <, <=, >, >=)
    for x in 0.0:0.1:1.0
        for y in 0.0:0.1:1.0
            @test f(x, y) ≈ discontinuity_eval(f, x, y)
        end
    end
end

@testset "approximation/majorization/minorization registrations" begin
    for f in (max, NaNMath.max)
        @test approximation_function(f) isa Function
        @test majorization_function(f) isa Function
    end
    for f in (min, NaNMath.min)
        @test approximation_function(f) isa Function
        @test minorization_function(f) isa Function
    end
    @test approximation_function(abs) isa Function
    @test majorization_function(abs) isa Function
    for f in (>=, >, <=, <, ==)
        @test approximation_function(f) isa Function
    end
end

@testset "majorization: $(nameof(f))" for (f, g) in (max => max, NaNMath.max => max)
    m = 100.0
    for a in -2.0:0.5:2.0, b in -2.0:0.5:2.0
        @test majorization_function(f)(m, a, b) >= g(a, b)
    end
end

@testset "minorization: $(nameof(f))" for (f, g) in (min => min, NaNMath.min => min)
    m = 100.0
    for a in -2.0:0.5:2.0, b in -2.0:0.5:2.0
        @test minorization_function(f)(m, a, b) <= g(a, b)
    end
end

@testset "majorization: abs" begin
    m = 100.0
    for a in -2.0:0.5:2.0
        @test majorization_function(abs)(m, a) >= abs(a)
    end
end

@testset "approximation convergence: $(nameof(f))" for (f, g) in (max => max, min => min, NaNMath.max => max, NaNMath.min => min)
    m = 1000.0
    for a in -2.0:0.5:2.0, b in -2.0:0.5:2.0
        isapprox(a, b; atol = 0.1) && continue
        @test approximation_function(f)(m, a, b) ≈ g(a, b) atol=1e-2
    end
end

@testset "approximation: comparison operators" begin
    m = 100.0
    for a in -2.0:0.5:2.0, b in -2.0:0.5:2.0
        # range [0, 1]
        @test 0.0 <= approximation_function(>=)(m, a, b) <= 1.0
        @test 0.0 <= approximation_function(>)(m, a, b) <= 1.0
        @test 0.0 <= approximation_function(<=)(m, a, b) <= 1.0
        @test 0.0 <= approximation_function(<)(m, a, b) <= 1.0
        # symmetry
        @test approximation_function(>=)(m, a, b) ≈ approximation_function(<=)(m, b, a)
        @test approximation_function(>)(m, a, b) ≈ approximation_function(<)(m, b, a)
    end
    # convergence
    m = 1000.0
    for a in -2.0:0.5:2.0, b in -2.0:0.5:2.0
        isapprox(a, b; atol = 0.1) && continue
        @test approximation_function(>=)(m, a, b) ≈ Float64(a >= b) atol=1e-2
        @test approximation_function(>)(m, a, b) ≈ Float64(a > b) atol=1e-2
        @test approximation_function(<=)(m, a, b) ≈ Float64(a <= b) atol=1e-2
        @test approximation_function(<)(m, a, b) ≈ Float64(a < b) atol=1e-2
    end
end

@testset "approximation: ==" begin
    m = 100.0
    for a in -2.0:0.5:2.0
        @test approximation_function(==)(m, a, a) ≈ 1.0 atol=1e-10
    end
    m = 1000.0
    for a in -2.0:0.5:2.0, b in -2.0:0.5:2.0
        abs(a - b) < 0.5 && continue
        @test approximation_function(==)(m, a, b) ≈ 0.0 atol=1e-2
    end
end
