using Symbolics
using SymbolicUtils, Test
using Symbolics: symtype, shape, wrap, unwrap, Unknown, Arr, arrterm, jacobian, @variables, value, get_variables, @arrayop
using Base: Slice
using SymbolicUtils: Sym, term, operation

@testset "arrays" begin
    @variables X[1:5, 1:5] Y[1:5, 1:5]
    @test_throws BoundsError X[1000]
    @test typeof(X) <: Arr
    @test shape(X) == Slice.((1:5, 1:5))
    @test shape(Y) == Slice.((1:5, 1:5))

    A = Y[2, :]
    @test typeof(A) <: Arr{Num,1}
    @test axes(A) == (1:5,)

    B = A[3:5]
    @test axes(B) == (Slice(1:3),)

    i = Sym{Int}(:i)
    j = Sym{Int}(:j)
    @test symtype(X[i, j]) == Real
    @test symtype(X[1, j]) == Real

    @variables t x[1:2](t)
    @test isequal(get_variables(0 ~ x[1]), [x[1]])
    @test Set(get_variables(2x)) == Set(collect(x)) # both array elements are present
    @test isequal(get_variables(2x[1]), [x[1]])
end

@testset "getindex" begin
    @variables X[1:5, 1:5] Y[1:5, 1:5]

    @test isequal(X[1, 1], wrap(term(getindex, unwrap(X), 1, 1)))

    XX = unwrap(X)
    @test isequal(unwrap(X[1, :]), Symbolics.@arrayop(XX[1, :], (j,), XX[1, j]))
    @test isequal(unwrap(X[:, 2]), Symbolics.@arrayop(XX[:, 2], (i,), XX[i, 2]))
    @test isequal(unwrap(X[:, 2:3]), Symbolics.@arrayop(XX[:, 2:3], (i, j), XX[i, j], (j in 2:3)))

    @variables t x[1:4](t)
    @syms i::Int
    @test isequal(x[i], operation(unwrap(x[i]))(t))
end

getdef(v) = getmetadata(v, Symbolics.VariableDefaultValue)
@testset "broadcast & scalarize" begin
    @variables A[1:5, 1:3] = 42 b[1:3] = [2, 3, 5] t x[1:4](t) u[1:1]
    AA = Symbolics.scalarize(A)
    bb = Symbolics.scalarize(b)
    @test all(isequal(42), getdef.(AA))
    @test getdef.(bb) == [2, 3, 5]
    @test isequal(Symbolics.scalarize([b .* 1; b .* 1]), [bb; bb])
    @test isequal(Symbolics.scalarize(b .^ 1), bb)
    c = A * b

    @test isequal(collect(sin.(x)),
        sin.([x[i] for i in 1:4]))

    @test isequal(Symbolics.scalarize(sin.(A .* c)[1, 1]),
        sin(A[1, 1] * (b[1] * A[1, 1] +
                       b[2] * A[1, 2] +
                       b[3] * A[1, 3])))

    D = Differential(t)
    @test isequal(collect(D.(x) ~ x), map(i -> D(x[i]) ~ x[i], eachindex(x)))
    @test_throws ArgumentError A ~ t

    # #448
    @test isequal(Symbolics.scalarize(u + u), [2u[1]])

    # #417
    @test isequal(Symbolics.scalarize(x', (1,1)), x[1])

    # #483
    # examples by @gronniger
    @variables A[1:2, 1:2]

    test_mat = [1 2; 3 4]
    repl_dict = Dict(Symbolics.scalarize(A .=> test_mat))
    A2 = A^2
    A3 = A^3
    A4 = A^4
    A5 = A^5
    A6 = A^6
    A7 = A^7

    @syms i::Int j::Int k::Int l::Int m::Int n::Int

    A_ = unwrap(A)
    A3_ = wrap(@arrayop (A_*A_*A_) (i, j) A_[i, k] * A_[k, l] * A_[l, j])
    A4_ = wrap(@arrayop (A_*A_*A_*A_) (i, j) A_[i, k] * A_[k, l] * A_[l, m] * A_[m, j])
    A5_ = wrap(@arrayop (A_*A_*A_*A_*A_) (i, j) A_[i, k] * A_[k, l] * A_[l, m] * A_[m, n] * A_[n, j])

    @test isequal(Symbolics.scalarize((A*A)[k,k]), A[k, 1]*A[1, k] + A[2, k]*A[k, 2])

    # basic tests:
    @test iszero((Symbolics.scalarize(A^2) * Symbolics.scalarize(A))[1,1] -
                  Symbolics.scalarize(A^3)[1,1])
    @testset "nested scalarize" begin
        @test isequal(substitute(Symbolics.scalarize(A2 ), repl_dict), test_mat^2)
        @test isequal(substitute(Symbolics.scalarize(A3_), repl_dict), test_mat^3)
        @test isequal(substitute(Symbolics.scalarize(A3 ), repl_dict), test_mat^3)
        @test isequal(substitute(Symbolics.scalarize(A4_), repl_dict), test_mat^4)
        @test isequal(substitute(Symbolics.scalarize(A4 ), repl_dict), test_mat^4)
        @test isequal(substitute(Symbolics.scalarize(A5_), repl_dict), test_mat^5)
        @test isequal(substitute(Symbolics.scalarize(A5 ), repl_dict), test_mat^5)
        @test isequal(substitute(Symbolics.scalarize(A6 ), repl_dict), test_mat^6)
        @test isequal(substitute(Symbolics.scalarize(A7 ), repl_dict), test_mat^7)
    end
    @test isequal(Symbolics.scalarize(x', (1, 1)), x[1])
end

@testset "Parent" begin
    @variables t x[1:4](t)
    x = unwrap(x)
    @test Symbolics.getparent(collect(x)[1]).metadata === x.metadata
end

n = 2
A = randn(n, n)
foo(x) = A * x # a function to represent symbolically, note, if this function is defined inside the testset, it's not found by the function fun_eval = eval(fun_ex)
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
        @syms foo(x::Symbolics.Arr{Num,1})::Symbolics.Arr{Num,1} # symbolic foo can not be created in global scope due to conflict with function foo
        foo(x) # return a symbolic call to foo
    end

    x0 = randn(n)
    @test foo(x0) == A * x0
    ex = symbolic_call(x)

    fun_genf = build_function(ex, x, expression=Val{false})
    @test_broken fun_genf(x0) == A * x0# UndefVarError: foo not defined

    # Generate an expression instead and eval it manually
    fun_ex = build_function(ex, x, expression=Val{true})
    fun_eval = eval(fun_ex)
    @test fun_eval(x0) == foo(x0)

    # Try to provide the hidden argument `expression_module` to solve the scoping issue
    @test_skip begin
        fun_genf = build_function(ex, x, expression=Val{false}, expression_module=Main) # UndefVarError: #_RGF_ModTag not defined
        fun_genf(x0) == A * x0
    end

    ## Jacobians
    @test Symbolics.value.(Symbolics.jacobian(foo(x), x)) == A
    @test_throws ErrorException Symbolics.value.(Symbolics.jacobian(ex, x))
end


@testset "Functions and Jacobians using manual @wrapped" begin
    @variables x[1:n]

    x0 = randn(n)
    @test foo(x0) == A * x0
    ex = foo(x)

    @test shape(ex) == shape(x)

    fun_genf = build_function(ex, x, expression=Val{false})
    @test fun_genf(x0) == A * x0

    # Generate an expression instead and eval it manually
    fun_ex = build_function(ex, x, expression=Val{true})
    fun_eval = eval(fun_ex)
    @test fun_eval(x0) == foo(x0)

    ## Jacobians
    @test value.(jacobian(foo(x), x)) == A
    @test value.(jacobian(ex, x)) == A
end

@testset "Rules" begin
    @variables X[1:10, 1:5] Y[1:5, 1:10] b[1:10]
    r = @rule ((~A * ~B) * ~C) => (~A * (~B * ~C)) where size(~A, 1) * size(~B, 2) >size(~B, 1)  * size(~C, 2)
    @test isequal(r(unwrap((X * Y) * b)), unwrap(X * (Y * b)))
end

@testset "2D Diffusion Composed With Stencil Interface" begin
    n = rand(8:32)

    @variables u[1:n, 1:n]
    @makearray v[1:n, 1:n] begin
        #interior
        v[2:end-1, 2:end-1] => @arrayop (i, j) u[i-1, j] + u[i+1, j] + u[i, j-1] + u[i, j+1] - 4 * u[i, j]
        #BCs
        v[1, 1:end] => 0.0
        v[n, 1:end] => 0.0
        v[1:end, 1] => 0.0
        v[1:end, n] => 0.0
    end

    #2D Diffusion composed
    @makearray ucx[1:n, 1:n] begin
        ucx[1:end, 1:end] => 0.0 # fill zeros
        ucx[2:end-1, 2:end-1] => @arrayop (i, j) u[i-1, j] + u[i+1, j] - 2 * u[i, j] (j in 2:n-1)
    end

    @makearray ucy[1:n, 1:n] begin
        ucy[1:end, 1:end] => 0.0 # fill zeros
        ucy[2:end-1, 2:end-1] => @arrayop (i, j) u[i, j-1] + u[i, j+1] - 2 * u[i, j] (i in 2:n-1)
    end

    uc = ucx .+ ucy

    global V, UC, UCX
    V, UC, UCX = v, uc, (ucx, ucy)
    @test isequal(collect(v), collect(uc))
end

@testset "ND Diffusion, Stencils with CartesianIndices" begin
    n = rand(8:32)
    N = 2

    @variables t u[fill(1:n, N)...](t)

    Igrid = CartesianIndices((fill(1:n, N)...))
    Iinterior = CartesianIndices((fill(2:n-1, N)...))

    function unitindices(N::Int) #create unit CartesianIndex for each dimension
        null = zeros(Int, N)
        if N == 0
            return CartesianIndex()
        else
            return map(1:N) do i
                unit_i = copy(null)
                unit_i[i] = 1
                CartesianIndex(Tuple(unit_i))
            end
        end
    end

    function Diffusion(N)
        ē = unitindices(N) # for i.e N = 3 => ē = [CartesianIndex((1,0,0)),CartesianIndex((0,1,0)),CartesianIndex((0,0,1))]

        Dss = map(1:N) do d
            @makearray u[Igrid] begin
                u[Iinterior] => @arrayop u[I-ē[d]] + u[I+ē[d]] - 2 * u[I]
            end
        end

        return reduce((D1, D2) -> D1 .+ D2, Dss)
    end

    D = Diffusion(N, n)

    Dxxu = @makearray u[1:n, 1:n] begin
        u[2:end-1, 2:end-1] => @arrayop (i, j) u[i-1, j] + u[i+1, j] - 2 * u[i, j]
    end

    Dyyu = @makearray u[1:n, 1:n] begin
        u[2:end-1, 2:end-1] => @arrayop (i, j) u[i, j-1] + u[i, j+1] - 2 * u[i, j]
    end

    @test D == Dxxu .+ Dyyu
end

@testset "Brusselator stencils" begin
    @variables t u[1:n, 1:n](t) v[1:n, 1:n](t)
    n = rand(8:32)

    limit(a, N) = a == N + 1 ? 1 : a == 0 ? N : a
    brusselator_f(x, y, t) = (((x - 0.3)^2 + (y - 0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.0

    x = y = range(0, stop=1, length=n)
    dx = step(x)

    A = 3.4
    alpha = 10.0

    @makearray dtu[1:n, 1:n] begin
        dtu[1:end, 1:end] => @arrayop (i, j) alpha * (u[limit(i - 1, n), j] + u[limit(i + 1, n), j] + u[i, limit(j + 1, n)] + u[i, limit(j - 1, n)] - 4u[i, j]) + 1.0 + u[i, j]^2 * v[i, j] - (A + 1) * u[i, j] + brusselator_f(x[i], y[j], t)
    end

    @makearray dtv[1:n, 1:n] begin
        dtv[1:end, 1:end] => @arrayop (i, j) alpha * (v[limit(i - 1, n), j] + v[limit(i + 1, n), j] + v[i, limit(j + 1, n)] + v[i, limit(j - 1, n)] - 4v[i, j]) - u[i, j]^2 * v[i, j] + A * u[i, j]
    end

    @makearray lapu[1:n, 1:n] begin
       lapu[1:end, 1:end] => @arrayop (i, j) u[i, j] + u[limit(i - 1, n), j] + u[limit(i + 1, n), j] + u[i, limit(j + 1, n)] + u[i, limit(j - 1, n)] - 4u[i, j]
    end

    @makearray lapv[1:n, 1:n] begin
        lapv[1:end, 1:end] => @arrayop (i, j) v[i, j] + v[limit(i - 1, n), j] + v[limit(i + 1, n), j] + v[i, limit(j + 1, n)] + v[i, limit(j - 1, n)] - 4v[i, j]
    end

    @makearray s[1:n, 1:n] begin
        s[1:end, 1:end] => @arrayop (i, j) brusselator_f(x[i], y[j], t)
    end

    @test dtu == 1 .+ v .* u.^2 .- (A + 1) .* u .+ alpha .* lapu .+ s
    @test dtv == A .* u .- u.^2 .* v .+ alpha .* lapv
end
