using Symbolics, Test
import SymbolicUtils: operation, arguments

signature_sum(a, b, c, d) = a + b + c + d
@register_symbolic signature_sum(a, b, c, d) true true begin
    signatures = [(Num, Num, Num, Num), (Num, Num, Real, Real)]
end

signature_empty(a, b) = a + b
@register_symbolic signature_empty(a, b) begin
    signatures = []
end

signature_array(a, b, c, d) = [a + b, c + d]
@register_array_symbolic signature_array(a, b, c, d) begin
    size = (2,)
    signatures = [(Num, Num, Num, Num), (Num, Num, Real, Real)]
end

const SIGNATURE_UNARY = Num
signature_unary(x) = sin(x)
@register_symbolic signature_unary(x) true begin
    signatures = [SIGNATURE_UNARY]
end

signature_wrapped(x::Real, y::Real) = x + y
Symbolics.@wrapped function signature_wrapped(x::Real, y::Real)
    x + y
end true begin
    signatures = [(Num, Real)]
end

signature_default(a, b) = a + b
@register_symbolic signature_default(a, b)

signature_array_default(a, b) = [a, b]
@register_array_symbolic signature_array_default(a, b) begin
    size = (2,)
end

signature_wrapped_default(x::Real) = x
Symbolics.@wrapped function signature_wrapped_default(x::Real)
    x
end

const SIGNATURE_PAIR = (Num, Real)
signature_pair(a, b) = a + b
@register_symbolic signature_pair(a, b) true true begin
    signatures = [SIGNATURE_PAIR]
end

signature_array_pair(a, b) = [a, b]
@register_array_symbolic signature_array_pair(a, b) begin
    size = (2,)
    signatures = [SIGNATURE_PAIR]
end

signature_wrapped_pair(x::Real, y::Real) = x + y
Symbolics.@wrapped function signature_wrapped_pair(x::Real, y::Real)
    x + y
end true begin
    signatures = [SIGNATURE_PAIR]
end

@testset "Explicit registration signatures" begin
    @variables a b c d
    @test length(methods(signature_pair)) == 2
    @test length(methods(signature_array_pair)) == 2
    @test length(methods(signature_wrapped_pair)) == 2
    @test operation(Symbolics.unwrap(signature_pair(a, b))) === signature_pair
    @test size(signature_array_pair(a, b)) == (2,)
    @test isequal(signature_wrapped_pair(a, b), a + b)
    @test length(methods(signature_default)) == 9
    @test length(methods(signature_array_default)) == 9
    @test length(methods(signature_wrapped_default)) == 3
    @test operation(Symbolics.unwrap(signature_default(a, b))) === signature_default
    @test size(signature_array_default(a, b)) == (2,)
    @test isequal(signature_wrapped_default(a), a)
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
    @test_throws ArgumentError Symbolics.registration_signatures(:([Num]))
    @test_throws ArgumentError Symbolics.registration_signatures(:(begin unknown = [] end))
    @test_throws ArgumentError Symbolics.registration_signatures(:(begin signatures end))
    body = :(validate_signature(x::Real, y::Real) = x + y)
    @test_throws ArgumentError Symbolics.wrap_func_expr(@__MODULE__, body, true, :([(Num, Num), (Num,)]))
    @test_throws ArgumentError Symbolics.wrap_func_expr(@__MODULE__, body, true, :([(Real, Real)]))
    @test_throws ArgumentError Symbolics.wrap_func_expr(@__MODULE__, body, true, :([(Num, 123)]))
    @test_throws ArgumentError Symbolics.wrap_func_expr(@__MODULE__, body, true, :([(Num, String)]))
end
