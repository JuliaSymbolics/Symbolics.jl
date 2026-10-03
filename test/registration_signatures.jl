using Symbolics, Test
import SymbolicUtils: operation, arguments

signature_sum(a, b, c, d) = a + b + c + d
@register_symbolic signature_sum(a, b, c, d) true true [(Num, Num, Num, Num), (Num, Num, Real, Real)]

signature_empty(a, b) = a + b
@register_symbolic signature_empty(a, b) true true []

signature_array(a, b, c, d) = [a + b, c + d]
@register_array_symbolic signature_array(a, b, c, d) begin
    size = (2,)
end true true [(Num, Num, Num, Num), (Num, Num, Real, Real)]

signature_unary(x) = sin(x)
@register_symbolic signature_unary(x) true true [Num]

signature_wrapped(x::Real, y::Real) = x + y
Symbolics.@wrapped function signature_wrapped(x::Real, y::Real)
    x + y
end true [(Num, Real)]

@testset "Explicit registration signatures" begin
    @variables a b c d
    @test length(methods(signature_sum)) == 3
    @test length(methods(signature_empty)) == 9
    @test signature_empty(1, 2) == 3
    @test operation(Symbolics.unwrap(signature_empty(a, b))) === signature_empty
    @test length(methods(signature_array)) == 3
    @test length(methods(signature_unary)) == 2
    @test length(methods(signature_wrapped)) == 2
    @test signature_sum(1, 2, 3, 4) == 10
    @test signature_array(1, 2, 3, 4) == [3, 7]
    @test signature_wrapped(1, 2) == 3
    for xs in ((a, b, c, d), (a, b, 3, 4), (a, b, c, 4))
        scalar = signature_sum(xs...)
        @test operation(Symbolics.unwrap(scalar)) === signature_sum
        @test all(x -> !(x isa Num), arguments(Symbolics.unwrap(scalar)))
        f = build_function(scalar, [a, b, c, d]; expression = Val(false))
        @test f([1, 2, 3, 4]) == 10
        array = signature_array(xs...)
        @test size(array) == (2,)
        f_array, _ = build_function(array, [a, b, c, d]; expression = Val(false))
        @test f_array([1, 2, 3, 4]) == [3, 7]
    end
    @test operation(Symbolics.unwrap(signature_unary(a))) === signature_unary
    @test isequal(signature_wrapped(a, b), a + b)
    body = :(validate_signature(x::Real, y::Real) = x + y)
    @test_throws ArgumentError Symbolics.wrap_func_expr(@__MODULE__, body, true, :([(Num, Num), (Num,)]))
    @test_throws ArgumentError Symbolics.wrap_func_expr(@__MODULE__, body, true, :([(Real, Real)]))
    @test_throws ArgumentError Symbolics.wrap_func_expr(@__MODULE__, body, true, :([(Num, 123)]))
    @test_throws ArgumentError Symbolics.wrap_func_expr(@__MODULE__, body, true, :([(Num, String)]))
end
