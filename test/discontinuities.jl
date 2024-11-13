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
