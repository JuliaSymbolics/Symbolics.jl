using Symbolics
using Symbolics: unwrap, operation, symtype
using Test

@variables x y

# A registered function that errors if it is ever actually evaluated. Used to detect whether
# the untaken branch of a conditional is being evaluated by generated code.
boom(v) = error("untaken branch was evaluated: boom($v)")
Symbolics.@register_symbolic boom(v)

@testset "Construction and operation identity" begin
    @test operation(unwrap(ifelse_eager(x > 0, x^2, 1 / x))) === ifelse_eager
    @test operation(unwrap(ifelse_branching(x > 0, x^2, 1 / x))) === ifelse_branching

    # condition must be Bool-typed, branches share type/shape with `ifelse`
    @test symtype(unwrap(ifelse_branching(x > 0, x, y))) ==
          symtype(unwrap(ifelse(x > 0, x, y)))

    # literal condition folds to the taken branch (matches `ifelse`)
    @test isequal(ifelse_branching(true, x, 1 / x), x)
    @test isequal(ifelse_eager(false, 1 / x, x), x)

    # literal (non-symbolic) branches are allowed
    @test operation(unwrap(ifelse_branching(x > 0, x, 5))) === ifelse_branching
    @test operation(unwrap(ifelse_eager(x > 0, 5, x))) === ifelse_eager
end

@testset "Both variants match `ifelse` numerically on valid inputs" begin
    for op in (ifelse_eager, ifelse_branching)
        f = build_function(op(x > 0, x^2, -x), x; expression = Val{false})
        @test f(3.0) == 9.0
        @test f(-3.0) == 3.0
    end
end

@testset "ifelse_branching does not evaluate the untaken branch" begin
    for cse in (true, false)
        f = build_function(
            ifelse_branching(x > 0, x^2 + 3x, boom(x)), x;
            expression = Val{false}, cse = cse
        )
        @test f(2.0) == 10.0          # takes the valid branch, `boom` never runs
        g = build_function(
            ifelse_branching(x > 0, boom(x), 1 / x), x;
            expression = Val{false}, cse = cse
        )
        @test g(-2.0) == -0.5
    end
end

# Counts how many times it is evaluated; used to detect duplicated computation.
const FIRE_COUNT = Ref(0)
fire(v) = (FIRE_COUNT[] += 1; v)
Symbolics.@register_symbolic fire(v)

@testset "multiply-referenced ifelse_branching is computed once under CSE" begin
    z = ifelse_branching(x > 0, fire(x)^2, boom(x))
    f = build_function(sin(z) + cos(z), x; expression = Val{false}, cse = true)
    FIRE_COUNT[] = 0
    @test f(2.0) == sin(4.0) + cos(4.0)
    @test FIRE_COUNT[] == 1     # the shared `if`/`else` ran once; `boom` not at all

    fvec = build_function([sin(z), cos(z)], x; expression = Val{false}, cse = true)[1]
    FIRE_COUNT[] = 0
    @test fvec(2.0) == [sin(4.0), cos(4.0)]
    @test FIRE_COUNT[] == 1
end

@testset "ifelse_eager evaluates both branches" begin
    for cse in (true, false)
        f = build_function(
            ifelse_eager(x > 0, x^2, boom(x)), x;
            expression = Val{false}, cse = cse
        )
        @test_throws ErrorException f(2.0)
    end
end

@testset "Nested and vector-valued conditionals stay lazy" begin
    f = build_function(
        ifelse_branching(x > 0, ifelse_branching(y > 0, x + y, boom(x)), boom(y)),
        x, y; expression = Val{false}, cse = true
    )
    @test f(2.0, 5.0) == 7.0

    fvec = build_function(
        [ifelse_branching(x > 0, x, boom(x)), 2y], x, y;
        expression = Val{false}, cse = true
    )[1]
    @test fvec(4.0, 3.0) == [4.0, 6.0]
end

@testset "Differentiation preserves the conditional variant" begin
    for op in (ifelse_eager, ifelse_branching)
        d = Symbolics.derivative(op(x > 0, x^2, x^3), x)
        @test operation(unwrap(d)) === op
        df = build_function(d, x; expression = Val{false})
        @test df(2.0) == 4.0        # d/dx x^2 = 2x
        @test df(-2.0) == 12.0      # d/dx x^3 = 3x^2
    end

    # the derivative of a branching conditional is itself lazy
    d = Symbolics.derivative(ifelse_branching(x > 0, x^2, boom(x)), x)
    df = build_function(d, x; expression = Val{false}, cse = true)
    @test df(2.0) == 4.0
end

@testset "substitute preserves the conditional variant" begin
    for op in (ifelse_eager, ifelse_branching)
        s = substitute(op(x > 0, x^2, 1 / x), Dict(x => y))
        @test operation(unwrap(s)) === op
    end
end

@testset "Lowering forms" begin
    @test Symbolics.toexpr(unwrap(ifelse_branching(x > 0, x, 1 / x))).head === :if
    @test Symbolics.toexpr(unwrap(ifelse_eager(x > 0, x, 1 / x))).head === :call

    # the eager call appears verbatim, the branching one as an `if`
    eager_code = build_function(ifelse_eager(x > 0, x, 1 / x), x; expression = Val{true})
    @test occursin("ifelse", string(eager_code))
    branch_code = string(build_function(
        ifelse_branching(x > 0, x, 1 / x), x; expression = Val{true}))
    @test occursin("if", branch_code)
end

@testset "Mixed scalar/array branches dispatch unambiguously" begin
    @variables v[1:3]
    for op in (ifelse_eager, ifelse_branching)
        # same-shape array branches build fine
        @test operation(unwrap(op(x > 0, v, v))) === op
        # mismatched scalar/array branches are shape-invalid: dispatch resolves to a
        # concrete method that then raises a shape error, rather than an ambiguity MethodError
        @test_throws ErrorException op(x > 0, x, v)
        @test_throws ErrorException op(x > 0, v, x)
    end
end

@testset "Non-Bool condition is rejected" begin
    # a non-Bool scalar condition with symbolic branches has no applicable method and errors
    @test_throws Exception ifelse_eager(1, x, y)
    @test_throws Exception ifelse_branching(1, x, y)
end

@testset "Sparsity detection matches `ifelse`" begin
    ref = Symbolics.hessian_sparsity(ifelse(x > 0, x^2 + y, x * y), [x, y])
    for op in (ifelse_eager, ifelse_branching)
        @test Symbolics.hessian_sparsity(op(x > 0, x^2 + y, x * y), [x, y]) == ref
    end
    refj = Symbolics.jacobian_sparsity([ifelse(x > 0, x * y, y)], [x, y])
    for op in (ifelse_eager, ifelse_branching)
        @test Symbolics.jacobian_sparsity([op(x > 0, x * y, y)], [x, y]) == refj
    end
end
