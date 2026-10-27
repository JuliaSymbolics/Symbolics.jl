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
