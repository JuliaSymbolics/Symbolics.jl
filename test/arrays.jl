using Symbolics
using SymbolicUtils, Test
using Symbolics: symtype, shape, wrap, unwrap, Unknown, Arr, arrterm, jacobian, @variables, value, get_variables
using Base: Slice
using SymbolicUtils: Sym, term, operation

@testset "arrays" begin
    @variables X[1:5, 1:5] Y[1:5, 1:5]
    @test_throws BoundsError X[1000]
    @test typeof(X) <: Arr
    @test shape(X) == Slice.((1:5, 1:5))
    @test shape(Y) == Slice.((1:5, 1:5))

    A = Y[2, :]
    @test typeof(A) <: Arr{Num, 1}
    @test axes(A) == (1:5,)

    B = A[3:5]
    @test axes(B) == (Slice(1:3),)

    i = Sym{Int}(:i)
    j = Sym{Int}(:j)
    @test symtype(X[i,j]) == Real
    @test symtype(X[1,j]) == Real

    @variables t x[1:2](t)
    @test isequal(get_variables(0 ~ x[1]), [x[1]])
    @test Set(get_variables(2x)) == Set(collect(x)) # both array elements are present
    @test isequal(get_variables(2x[1]), [x[1]])
end

@testset "getindex" begin
    @variables X[1:5, 1:5] Y[1:5, 1:5]

    @test isequal(X[1,1], wrap(term(getindex, unwrap(X), 1,1)))

    XX = unwrap(X)
    @test isequal(unwrap(X[1, :]), Symbolics.@arrayop(XX[1,:], (j,), XX[1, j]))
    @test isequal(unwrap(X[:, 2]), Symbolics.@arrayop(XX[:,2], (i,), XX[i, 2]))
    @test isequal(unwrap(X[:, 2:3]), Symbolics.@arrayop(XX[:,2:3], (i,j), XX[i, j], (+), (j in 2:3)))

    @variables t x[1:4](t)
    @syms i::Int
    @test isequal(x[i], operation(unwrap(x[i]))(t))
end

getdef(v) = getmetadata(v, Symbolics.VariableDefaultValue)
@testset "broadcast & scalarize" begin
    @variables A[1:5,1:3]=42 B[1:2, 1:2] b[1:3]=[2, 3, 5] t x[1:4](t) u[1:1]
    AA = Symbolics.scalarize(A)
    bb = Symbolics.scalarize(b)
    @test all(isequal(42), getdef.(AA))
    @test getdef.(bb) == [2, 3, 5]
    @test isequal(Symbolics.scalarize([b.*1; b.*1]), [bb; bb])
    @test isequal(Symbolics.scalarize(b.^1), bb)
    c = A*b

    @test isequal(collect(sin.(x)),
                  sin.([x[i] for i in 1:4]))

    @test isequal(Symbolics.scalarize(sin.(A.*c)[1,1]),
                  sin(A[1, 1]*(b[1]*A[1, 1] +
                               b[2]*A[1, 2] +
                               b[3]*A[1, 3])))

    D = Differential(t)
    @test isequal(collect(D.(x) ~ x), map(i->D(x[i]) ~ x[i], eachindex(x)))
    @test_throws ArgumentError A ~ t

    # #448
    @test isequal(Symbolics.scalarize(u + u), [2u[1]])

    # #417
    @test isequal(Symbolics.scalarize(x', (1,1)), x[1])

    test_mat = [1 2; 3 4]
    repl_dict = Dict(Symbolics.scalarize(B .=> test_mat))
    @test isequal(substitute(Symbolics.scalarize(B^2), repl_dict), test_mat^2)
    @test isequal(substitute(Symbolics.scalarize(B^3), repl_dict), test_mat^3)
    @test isequal(substitute(Symbolics.scalarize(B^4), repl_dict), test_mat^4)
    @test_broken isequal(substitute(Symbolics.scalarize(B^5), repl_dict), test_mat^5)
    @test_broken isequal(substitute(Symbolics.scalarize(B^6), repl_dict), test_mat^6)
    @test_broken isequal(substitute(Symbolics.scalarize(B^7), repl_dict), test_mat^7)
end

@testset "Parent" begin
    @variables t x[1:4](t)
    x = unwrap(x)
    @test Symbolics.getparent(collect(x)[1]).metadata === x.metadata
end

n = 2
A = randn(n,n)
foo(x) = A*x # a function to represent symbolically, note, if this function is defined inside the testset, it's not found by the function fun_eval = eval(fun_ex)
function Symbolics.propagate_ndims(::typeof(foo), x)
    ndims(x)
end
function Symbolics.propagate_shape(::typeof(foo), x)
    shape(x)
end
@wrapped function foo(x::AbstractVector)
    t = arrterm(foo, x)
    setmetadata(t, Symbolics.ScalarizeCache, Ref{Any}(nothing))
end

#=
The following two testsets test jacobians for symbolic functions of symbolic arrays. Getting it to work currently requires the user to manually specify propagate_ndims and propagate_shape. Making use of `@register` or `@syms` to specify symbolic functions fail to propagate array shape
=#

@testset "Functions and Jacobians using @syms" begin
    @variables x[1:n]

    function symbolic_call(x)
        @syms foo(x::Symbolics.Arr{Num, 1})::Symbolics.Arr{Num, 1} # symbolic foo can not be created in global scope due to conflict with function foo
        foo(x) # return a symbolic call to foo
    end

    x0 = randn(n)
    @test foo(x0) == A*x0
    ex = symbolic_call(x)

    fun_genf = build_function(ex, x, expression=Val{false})
    @test_broken fun_genf(x0) == A*x0# UndefVarError: foo not defined

    # Generate an expression instead and eval it manually
    fun_ex = build_function(ex, x, expression=Val{true})
    fun_eval = eval(fun_ex)
    @test fun_eval(x0) == foo(x0)

    # Try to provide the hidden argument `expression_module` to solve the scoping issue
    @test_skip begin
        fun_genf = build_function(ex, x, expression=Val{false}, expression_module=Main) # UndefVarError: #_RGF_ModTag not defined
        fun_genf(x0) == A*x0
    end

    ## Jacobians
    @test Symbolics.value.(Symbolics.jacobian(foo(x), x)) == A
    @test_skip Symbolics.value.(Symbolics.jacobian(ex , x)) == A #ERROR: axes of foo(x[1:2]) not known
end


@testset "Functions and Jacobians using manual @wrapped" begin
    @variables x[1:n]

    x0 = randn(n)
    @test foo(x0) == A*x0
    ex = foo(x)

    @test shape(ex) == shape(x)

    fun_genf = build_function(ex, x, expression=Val{false})
    @test fun_genf(x0) == A*x0

    # Generate an expression instead and eval it manually
    fun_ex = build_function(ex, x, expression=Val{true})
    fun_eval = eval(fun_ex)
    @test fun_eval(x0) == foo(x0)

    ## Jacobians
    @test value.(jacobian(foo(x), x)) == A
    @test value.(jacobian(ex , x)) == A
end
