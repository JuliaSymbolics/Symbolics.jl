using Symbolics, SparseArrays, Test
using Symbolics: unwrap, is_derivative

struct OpaqueOp
    A::SparseMatrixCSC{Float64, Int}
end
apply_op(A::OpaqueOp, x::AbstractVector) = A.A * x
@register_array_symbolic apply_op(A::OpaqueOp, x::AbstractVector) begin
    size = size(x)
    eltype = eltype(x)
    ndims = 1
end

offset(A::AbstractMatrix, x) = x + 1
@register_symbolic offset(A::AbstractMatrix, x)
@register_derivative offset(A, x) 2 (x - x + 1)

N = 4
@variables u[1:N]
@variables Kp::OpaqueOp
@variables Mp::SparseMatrixCSC{Float64, Int}

@testset "array-valued call without derivative rule" begin
    ex = Symbolics.scalarize(apply_op(Kp, u))
    @test nnz(Symbolics.jacobian_sparsity(ex, collect(u))) == N^2
    J = Symbolics.jacobian(ex, collect(u))
    # unknown partials stay as unexpanded derivatives, never silent zeros
    @test !any(iszero, J)
    @test all(j -> is_derivative(unwrap(j)), J)
    Jsp = Symbolics.sparsejacobian(ex, collect(u))
    @test nnz(Jsp) == N^2
    @test_throws Symbolics.DerivativeNotDefinedError expand_derivatives(
        Differential(u[1])(ex[1]); throw_no_derivative = true
    )
end

@testset "array-typed symbol used as a scalar coefficient" begin
    ex = Symbolics.scalarize(Mp * u)
    @test_throws ArgumentError Symbolics.jacobian_sparsity(ex, collect(u))
    @test_throws ArgumentError Symbolics.jacobian(ex, collect(u))
    @test_throws ArgumentError Symbolics.sparsejacobian(ex, collect(u))
    # independent of the variables: no error
    @variables p
    @test isequal(Symbolics.jacobian([Mp * p + u[1]], collect(u))[1, 1], 1)
    # operand of a registered scalar function with a derivative rule: fine
    J = Symbolics.jacobian([offset(Mp, u[1])], collect(u))
    @test isequal(J[1, 1], 1) && all(iszero, J[1, 2:end])
end

@register_array_symbolic foo(x) begin
    size = (2,)
    eltype = Real
    ndims = 1
end
foo(x) = [x, 2x]

@testset "explicit element dependencies and unrelated elements" begin
    @variables v[1:2]
    J = Symbolics.jacobian([foo(v[1])[1]], collect(v))
    @test !iszero(J[1, 1]) && iszero(J[1, 2])
    @test Symbolics.jacobian_sparsity([foo(v[1])[1]], collect(v)) == sparse([1], [1], true, 1, 2)
    # `Mp * u[2]` does not depend on `u[1]`
    @test iszero(Symbolics.jacobian([Mp * u[2]], [u[1]])[1, 1])
    @test nnz(Symbolics.jacobian_sparsity([Mp * u[2]], [u[1]])) == 0
    @test_throws ArgumentError Symbolics.jacobian([Mp * u[2]], [u[2]])
end

@register_array_symbolic bar(x::AbstractVector) begin
    size = (2,)
    eltype = Real
    ndims = 1
end
bar(x::AbstractVector) = x

@testset "nested opaque calls, element arguments and call-valued arrays" begin
    @variables v[1:2] t w(t)[1:2]
    nested = bar(foo(v[1]))[1]
    @test !iszero(Symbolics.jacobian([nested], collect(v))[1, 1])
    @test iszero(Symbolics.jacobian([nested], collect(v))[1, 2])
    @test is_derivative(unwrap(expand_derivatives(Differential(v[1])(nested))))
    @test iszero(expand_derivatives(Differential(v[2])(nested)))
    @test Symbolics.jacobian_sparsity([nested], collect(v)) == sparse([1], [1], true, 1, 2)
    @test !iszero(Symbolics.jacobian([bar(v)[1]], collect(v))[1, 2])
    # distinct elements of u(t)[1:2] are independent
    @test iszero(Symbolics.jacobian([Mp * w[2]], [w[1]])[1, 1])
    @test nnz(Symbolics.jacobian_sparsity([Mp * w[2]], [w[1]])) == 0
    @test_throws ArgumentError Symbolics.jacobian([Mp * w[2]], [w[2]])
    @test_throws ArgumentError Symbolics.jacobian_sparsity([Mp * w[1]], collect(w))
end

ignore_argument(A, x) = x + 1
@register_symbolic ignore_argument(A, x)
@register_derivative ignore_argument(A, x) 1 (x - x)
@register_derivative ignore_argument(A, x) 2 (x - x + 1)

use_argument(A, x) = x + 1
@register_symbolic use_argument(A, x)
@register_derivative use_argument(A, x) 1 x
@register_derivative use_argument(A, x) 2 (x - x + 1)

@testset "registered zero partial skips the opaque operand check" begin
    ex = [ignore_argument(Mp * u[1], u[2])]
    J = Symbolics.jacobian(ex, collect(u))
    @test iszero(J[1, 1]) && isequal(J[1, 2], 1)
    Jsp = Symbolics.sparsejacobian(ex, collect(u))
    @test Jsp[1, 2] == 1 && iszero(Jsp[1, 1])
    @test Symbolics.jacobian_sparsity(ex, collect(u))[1, 2]
    # negative control: a nonzero partial wrt the arithmetic-array operand still throws
    ex = [use_argument(Mp * u[1], u[2])]
    @test_throws ArgumentError Symbolics.jacobian(ex, collect(u))
    @test_throws ArgumentError Symbolics.sparsejacobian(ex, collect(u))
    @test_throws ArgumentError Symbolics.jacobian_sparsity(ex, collect(u))
end

err_partial(z) = z
@register_symbolic err_partial(z)
@register_derivative err_partial(z) 1 error("unsupported partial")

err_partial_arr(A, x) = x
@register_symbolic err_partial_arr(A, x)
@register_derivative err_partial_arr(A, x) 1 error("unsupported partial")
@register_derivative err_partial_arr(A, x) 2 (x - x + 1)

@testset "registered partials are only evaluated when an opaque operand needs checking" begin
    @variables y
    # no opaque operand: the erroring rule is never evaluated
    @test iszero(Symbolics.jacobian([err_partial(y)], [u[1]])[1, 1])
    @test iszero(Symbolics.sparsejacobian([err_partial(y)], [u[1]])[1, 1])
    @test iszero(nnz(Symbolics.jacobian_sparsity([err_partial(y)], [u[1]])))
    @test nnz(Symbolics.jacobian_sparsity([err_partial(u[1])], [u[1]])) == 1
    # opaque operand and an erroring rule: the original ArgumentError, not the rule's error
    ex = [err_partial_arr(Mp * u[1], u[2])]
    @test_throws ArgumentError Symbolics.jacobian(ex, collect(u))
    @test_throws ArgumentError Symbolics.sparsejacobian(ex, collect(u))
    @test_throws ArgumentError Symbolics.jacobian_sparsity(ex, collect(u))
end

ignore_counter(A, x) = x + 1
@register_symbolic ignore_counter(A, x)
@register_derivative ignore_counter(A, x) 1 (x - x)
@register_derivative ignore_counter(A, x) 2 (x - x + 1)

sum_counter(A, x) = sum(A) + x
@register_symbolic sum_counter(A, x)
@register_derivative sum_counter(A, x) 2 (x - x + 1)

@variables Mq::SparseMatrixCSC{Float64, Int} p

@testset "zero-partial calls nested in every direction around opaque arithmetic" begin
    function check(expr)
        J = Symbolics.jacobian([expr], collect(u))
        @test iszero(J[1, 1]) && isequal(J[1, 2], 1)
        Jsp = Symbolics.sparsejacobian([expr], collect(u))
        @test iszero(Jsp[1, 1]) && Jsp[1, 2] == 1
        @test Symbolics.jacobian_sparsity([expr], collect(u))[1, 2]
    end
    # zero-partial call inside the opaque arithmetic
    check(sum_counter(Mp * ignore_counter(u[1], p), u[2]))
    # several levels of nesting
    check(sum_counter(Mp * ignore_counter(Mq * ignore_counter(u[1], p), p), u[2]))
    # zero-partial calls both inside and enclosing the opaque arithmetic
    check(ignore_argument(Mp * ignore_counter(u[1], p), u[2]))

    # negative controls: still rejected
    for expr in (
        sum_counter(Mp * use_argument(u[1], p), u[2]),
        sum_counter(Mp * err_partial_arr(u[1], p), u[2]),
        sum_counter(Mp * ignore_counter(u[2], u[1]), u[2]),
        Mp * ignore_counter(u[1], p) + Mp * u[1],
    )
        @test_throws ArgumentError Symbolics.jacobian([expr], collect(u))
        @test_throws ArgumentError Symbolics.sparsejacobian([expr], collect(u))
        @test_throws ArgumentError Symbolics.jacobian_sparsity([expr], collect(u))
    end
end
