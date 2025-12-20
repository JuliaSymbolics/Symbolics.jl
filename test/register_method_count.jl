using Test
using Symbolics

"""
Count the number of function definitions in a macro-expanded expression.
"""
function count_function_defs(expr)
    count = Ref(0)
    _count_function_defs!(expr, count)
    return count[]
end

function _count_function_defs!(expr, count)
    if expr isa Expr
        if expr.head === :function
            count[] += 1
        end
        for arg in expr.args
            _count_function_defs!(arg, count)
        end
    end
end

@testset "Method count scaling for @register_symbolic" begin
    # Test that the number of methods generated scales only with UNTYPED arguments.
    #
    # Structure of generated code:
    # - 1 impl function (the actual implementation)
    # - N dispatch functions (where N = prod(type_options per arg) - 1)
    # - 1 promote_shape function
    #
    # Type options per argument:
    # - Typed args (e.g., ::Float64) get 1 type option (the specified type)
    # - Untyped args (defaulting to ::Real) get 3 type options: (Real, BasicSymbolic, Num)
    #
    # For k typed + m untyped args: 1^k * 3^m - 1 dispatch methods

    # Case 1: 2 untyped args
    # Type options: 3 × 3 = 9 combinations, minus (Real, Real) = 8 dispatch
    # Total: 8 dispatch + 1 impl + 1 promote_shape = 10
    f2u(a, b) = a * b
    expr2u = @macroexpand @register_symbolic f2u(a, b)
    n2u = count_function_defs(expr2u)
    @test n2u == 10

    # Case 2: 1 typed + 1 untyped
    # Type options: 1 × 3 - 1 = 2 dispatch
    # Total: 2 dispatch + 1 impl + 1 promote_shape = 4
    f1t1u(a::Float64, b) = a * b
    expr1t1u = @macroexpand @register_symbolic f1t1u(a::Float64, b)
    n1t1u = count_function_defs(expr1t1u)
    @test n1t1u == 4

    # Case 3: 2 typed + 2 untyped
    # Type options: 1 × 1 × 3 × 3 - 1 = 8 dispatch
    # Total: 8 dispatch + 1 impl + 1 promote_shape = 10
    f2t2u(a::Float64, b::Float64, c, d) = a * b * c * d
    expr2t2u = @macroexpand @register_symbolic f2t2u(a::Float64, b::Float64, c, d)
    n2t2u = count_function_defs(expr2t2u)
    @test n2t2u == 10

    # Case 4: 6 typed + 2 untyped
    # Type options: 1^6 × 3^2 - 1 = 8 dispatch
    # Total: 8 dispatch + 1 impl + 1 promote_shape = 10
    f6t2u(a::Float64, b::Float64, c::Float64, d::Float64, e::Float64, f::Float64, g, h) =
        a * b * c * d * e * f * g * h
    expr6t2u = @macroexpand @register_symbolic f6t2u(
        a::Float64, b::Float64, c::Float64, d::Float64, e::Float64, f::Float64, g, h)
    n6t2u = count_function_defs(expr6t2u)
    @test n6t2u == 10

    # Case 5: All typed (4 args)
    # Type options: 1^4 - 1 = 0 dispatch
    # Total: 0 dispatch + 1 impl + 1 promote_shape = 2
    f4t(a::Float64, b::Float64, c::Float64, d::Float64) = a * b * c * d
    expr4t = @macroexpand @register_symbolic f4t(
        a::Float64, b::Float64, c::Float64, d::Float64)
    n4t = count_function_defs(expr4t)
    @test n4t == 2
end
