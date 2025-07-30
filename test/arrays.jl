using Symbolics
using SymbolicUtils, Test
using Symbolics: symtype, shape, wrap, unwrap, Unknown, Arr, array_term, jacobian, @variables, value, get_variables, @arrayop, getname, metadata, scalarize
using Base: Slice
using SymbolicUtils: Sym, term, operation
import LinearAlgebra: dot
import ..limit2

struct TestMetaT end
Symbolics.option_to_metadata_type(::Val{:test_meta}) = TestMetaT

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

    @variables t x(t)[1:2]
    @test isequal(get_variables(0 ~ x[1]), [x[1]])
    @test Set(get_variables(2x)) == Set(collect(x)) # both array elements are present
    @test isequal(get_variables(2x[1]), [x[1]])
end

@testset "getname" begin
    @variables t x(t)[1:4]
    v = Symbolics.lower_varname(unwrap(x[2]), unwrap(t), 2)
    @test operation(v) == getindex
    @test arguments(v)[2] == 2
    @test getname(v) == getname(arguments(v)[1]) == Symbol("xˍtt")
end

@testset "getindex" begin
    @variables X[1:5, 1:5] Y[1:5, 1:5]

    @test isequal(X[1, 1], wrap(term(getindex, unwrap(X), 1, 1)))

    XX = unwrap(X)
    @test isequal(unwrap(X[1, :]), Symbolics.@arrayop((j,), XX[1, j], term=XX[1, :]))
    @test isequal(unwrap(X[:, 2]), Symbolics.@arrayop((i,), XX[i, 2], term=XX[:, 2]))
    @test isequal(unwrap(X[:, 2:3]), Symbolics.@arrayop((i, j), XX[i, j], (j in 2:3), term=XX[:, 2:3]))

    @variables t x(t)[1:4]
    @syms i::Int
    @test isequal(x[i], operation(unwrap(x))(t)[i])

    # https://github.com/JuliaSymbolics/Symbolics.jl/issues/842
    # getindex should keep metadata
    @variables tv v(tv)[1:2] [test_meta = 4] v2(tv)[1:3] [test_meta=[1, 2, 3]]
    @test !isnothing(metadata(unwrap(v)))
    @test hasmetadata(unwrap(v), TestMetaT)
    @test getmetadata(unwrap(v), TestMetaT) == 4
    @test getmetadata(unwrap(v2), TestMetaT) == [1, 2, 3]
    vs = scalarize(v)
    vsw = unwrap.(vs)
    vs2 = scalarize(v2)
    vsw2 = unwrap.(vs2)
    @test !isnothing(metadata(vsw[1]))
    @test hasmetadata(vsw[1], TestMetaT)
    @test getmetadata(vsw[1], TestMetaT) == 4
    @test getmetadata.(vsw2, TestMetaT) == [1, 2, 3]
    @test !isnothing(metadata(unwrap(v[1])))
    @test hasmetadata(unwrap(v[1]), TestMetaT)
    @test getmetadata(unwrap(v[1]), TestMetaT) == 4
end

@testset "maketerm" begin
    @variables A[1:5, 1:5] B[1:5, 1:5] C

    T = unwrap(3A)
    @test isequal(T, Symbolics.maketerm(typeof(T), operation(T), arguments(T), nothing))
    T2 = unwrap(3B)
    @test isequal(T2, Symbolics.maketerm(typeof(T), operation(T), [*, 3, unwrap(B)], nothing))
    T3 = unwrap(A .^ 2)
    @test isequal(T3, Symbolics.maketerm(typeof(T3), operation(T3), arguments(T3), nothing))
    T4 = unwrap(A .* C)
    @test isequal(T4, Symbolics.maketerm(typeof(T4), operation(T4), arguments(T4), nothing))
end

getdef(v) = getmetadata(v, Symbolics.VariableDefaultValue)
@testset "broadcast & scalarize" begin
    @variables A[1:5,1:3]=42 b[1:3]=[2, 3, 5] t x(t)[1:4] u[1:1]
    AA = Symbolics.scalarize(A)
    bb = Symbolics.scalarize(b)
    @test all(isequal(42), getdef.(AA))
    @test getdef.(bb) == [2, 3, 5]
    @test isequal(Symbolics.scalarize([b .* 1; b .* 1]), [bb; bb])
    @test isequal(Symbolics.scalarize(b .^ 1), bb)
    c = A * b

    # Test for issue #575: adjoint multiplication with Symbolics.Arr should not be ambiguous
    @variables d[1:3] E[1:3, 1:3]
    d_vec = collect(d)  # Convert to Vector{Num}
    # These should work without MethodError due to ambiguity
    @test d' isa Adjoint{Num, Vector{Num}}
    @test E isa Symbolics.Arr{Num, 2}
    result1 = d' * E  # This was causing ambiguity error
    result2 = d' * inv(E) * d  # The original failing expression from issue #575
    @test size(result1) == (1, 3)
    @test size(result2) == (1, 1)

    @test isequal(collect(sin.(x)),
        sin.([x[i] for i in 1:4]))

    @test isequal(Symbolics.scalarize(sin.(A .* c)[1, 1]),
        sin(A[1, 1] * (b[1] * A[1, 1] +
                       b[2] * A[1, 2] +
                       b[3] * A[1, 3])))

    D = Differential(t)
    @test isequal(collect(D.(x) .~ x), map(i -> D(x[i]) ~ x[i], eachindex(x)))
    @test_throws ArgumentError A ~ t
    @test isequal(D(x[1]), D(x)[1])
    a = Symbolics.unwrap(D(x)[1])
    @test Symbolics.operation(a) == D
    @test isequal(only(Symbolics.arguments(a)), Symbolics.unwrap(x[1]))

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
    A3_ = wrap(@arrayop (i, j) A_[i, k] * A_[k, l] * A_[l, j])
    A4_ = wrap(@arrayop (i, j) A_[i, k] * A_[k, l] * A_[l, m] * A_[m, j])
    A5_ = wrap(@arrayop (i, j) A_[i, k] * A_[k, l] * A_[l, m] * A_[m, n] * A_[n, j])

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

    ##653
    Symbolics.scalarize(inv(A)[1,1])

    ##895
    @test inv(Num.(reshape([1],1,1)))==Num.(reshape([1],1,1))

    # #831
    @syms symT sym1(symT) sym2(symT)
    symvec = [sym1(symT), sym2(symT)]
    weights = rand(2)
    @test isequal(symvec'weights, dot(symvec, weights))
    @test isequal(weights'symvec, dot(weights, symvec))
    @syms sym3(symT)::Real sym4(symT)::Real
    symvec = [sym3(symT), sym4(symT)]
    @test isequal(symvec'weights, weights'symvec)

    # ModelingToolkit.jl#1736
    #
    @variables t F(t)[1:1]

    @test isequal(collect(F ./ t), [F[1] / t])
end

n = 2
A = randn(n, n)
foo(x) = A * x # a function to represent symbolically, note, if this function is defined inside the testset, it's not found by the function fun_eval = eval(fun_ex)
@register_array_symbolic foo(x::Vector{Real}) begin
    size = (n,)
    eltype = eltype(x)
    ndims = 1
end

#=
The following two testsets test jacobians for symbolic functions of symbolic arrays. Getting it to work currently requires the user to manually specify propagate_ndims and propagate_shape. Making use of `@register` or `@syms` to specify symbolic functions fail to propagate array shape
=#

@testset "Functions and Jacobians using @syms" begin
    @variables x[1:n]

    x0 = randn(n)
    @test foo(x0) == A * x0
    ex = foo(x)

    fun_oop, fun_iip = build_function(ex, x, expression=Val{false})
    @test fun_oop(x0) == A * x0# UndefVarError: foo not defined

    # Generate an expression instead and eval it manually
    fun_ex_oop, fun_ex_iip = build_function(ex, x, expression=Val{true})
    fun_eval = eval(fun_ex_oop)
    @test fun_eval(x0) == foo(x0)
end


@testset "Functions and Jacobians using manual @wrapped" begin
    @variables x[1:n]

    x0 = randn(n)
    @test foo(x0) == A * x0
    ex = foo(x)

    @test shape(ex) == shape(x)

    fun_oop, fun_iip = build_function(ex, x, expression=Val{false})
    @test fun_oop(x0) == A * x0

    # Generate an expression instead and eval it manually
    fun_ex_oop, fun_ex_ip = build_function(ex, x, expression=Val{true})
    fun_eval = eval(fun_ex_oop)
    @test fun_eval(x0) == foo(x0)

    ## Jacobians
    @test_broken value.(jacobian(foo(x), x)) == A
    @test_broken value.(jacobian(ex , x)) == A

end

@testset "Rules" begin
    @variables X[1:10, 1:5] Y[1:5, 1:10] b[1:10]
   #r = @rule ((~A * ~B) * ~C) => (~A * (~B * ~C)) where (size(~A, 1) * size(~B, 2) >size(~B, 1)  * size(~C, 2))
   #@test isequal(r(unwrap((X * Y) * b)), unwrap(X * (Y * b)))
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

    @variables t u(t)[fill(1:n, N)...]

    Igrid = CartesianIndices((fill(1:n, N)...,))
    Iinterior = CartesianIndices((fill(2:n-1, N)...,))

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
    function Diffusion(N, n)
        ē = unitindices(N) # for i.e N = 3 => ē = [CartesianIndex((1,0,0)),CartesianIndex((0,1,0)),CartesianIndex((0,0,1))]

        Dss = map(1:N) do d
            ranges = CartesianIndices((map(i->d == i ? (2:n-1) : (1:n), 1:N)...,))
            @makearray x[1:n, 1:n] begin
                x[1:n, 1:n] => 0
                x[ranges] => @arrayop (i, j) u[CartesianIndex(i, j)-ē[d]] +
                                                u[CartesianIndex(i, j)+ē[d]] - 2 * u[i, j]
            end
        end

        return reduce((D1, D2) -> D1 .+ D2, Dss)
    end

    D = Diffusion(N, n)

    @makearray Dxxu[1:n, 1:n] begin
        Dxxu[1:n, 1:n] => 0
        Dxxu[2:end-1, 1:end] => @arrayop (i, j) u[i-1, j] + u[i+1, j] - 2 * u[i, j]
    end

    @makearray Dyyu[1:n, 1:n] begin
        Dyyu[1:n, 1:n] => 0
        Dyyu[1:end, 2:end-1] => @arrayop (i, j) u[i, j-1] + u[i, j+1] - 2 * u[i, j]
    end

    @test isequal(collect(D), collect(Dxxu .+ Dyyu))
end

@testset "Brusselator stencil" begin
    n = 8
    @variables t u(t)[1:n, 1:n] v(t)[1:n, 1:n]

    brusselator_f(x, y, t) = (((x - 0.3)^2 + (y - 0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.0

    x = y = range(0, stop=1, length=n)
    dx = step(x)

    A = 3.4
    alpha = 10.0

    dtu = @arrayop (i, j) alpha * (u[limit2(i - 1, n), j] +
                                   u[limit2(i + 1, n), j] +
                                   u[i, limit2(j + 1, n)] +
                                   u[i, limit2(j - 1, n)] -
                                   4u[i, j]) +
                          1.0 + u[i, j]^2 * v[i, j] - (A + 1) *
                            u[i, j] + brusselator_f(x[i], y[j], t) i in 1:n j in 1:n
    dtv = @arrayop (i, j) alpha * (v[limit2(i - 1, n), j] +
                                   v[limit2(i + 1, n), j] +
                                   v[i, limit2(j + 1, n)] +
                                   v[i, limit2(j - 1, n)] -
                                   4v[i, j]) -
                          u[i, j]^2 * v[i, j] + A * u[i, j] i in 1:n j in 1:n
    lapu = @arrayop (i, j) (u[limit2(i - 1, n), j] +
                            u[limit2(i + 1, n), j] +
                            u[i, limit2(j + 1, n)] +
                            u[i, limit2(j - 1, n)] -
                            4u[i, j]) i in 1:n j in 1:n
    lapv = @arrayop (i, j) (v[limit2(i - 1, n), j] +
                            v[limit2(i + 1, n), j] +
                            v[i, limit2(j + 1, n)] +
                            v[i, limit2(j - 1, n)] -
                            4v[i, j]) i in 1:n j in 1:n
    s = brusselator_f.(x, y', t)

    dtu = wrap(dtu)
    dtv = wrap(dtv)
    lapu = wrap(lapu)
    lapv = wrap(lapv)

    g, f = build_function(dtu, u, v, t, expression=Val{false}, nanmath = false)
    du = zeros(Num, 8, 8)
    f(du, u,v,t)
    @test isequal(collect(du), collect(dtu))

    @test isequal(collect(dtu), collect(1 .+ v .* u.^2 .- (A + 1) .* u .+ alpha .* lapu .+ s))
    @test isequal(collect(dtv), collect(A .* u .- u.^2 .* v .+ alpha .* lapv))
end

@testset "Unwrapped array equality" begin
    @variables x[1:3]
    ux = unwrap(x)
    @test isequal(x, x)
    @test isequal(x, ux)
    @test isequal(ux, x)
end

@testset "Array expression substitution" begin
    @variables x[1:3] p[1:3, 1:3]
    bar(x, p) = p * x
    @register_array_symbolic bar(x::AbstractVector, p::AbstractMatrix) begin
        size = size(x)
        eltype = promote_type(eltype(x), eltype(p))
        ndims = 1
    end

    @test isequal(substitute(bar(x, p), x => ones(3)), bar(ones(3), p))
    @test isequal(substitute(bar(x, p), Dict(x => ones(3), p => ones(3, 3))), wrap(3ones(3)))
    @test isequal(substitute(bar(x, p), [x => ones(3), p => ones(3, 3)]), wrap(3ones(3)))
end

@testset "Partial array substitution" begin
    @variables x[1:3] A[1:2, 1:2, 1:2]

    @test substitute(x[1], Dict(x => [1, 2, 3])) === Num(1)
    @test substitute(A[1,2,1], Dict(A => reshape(1:8, 2, 2, 2))) === Num(3)

    @test substitute(A[1,2,1], Dict(A[1,2,1] => 9)) === Num(9)
end

@testset "Hashes" begin
    @variables u[1:7]
    a, b, c = u[1:5], u[2:6], u[3:7]
    @test !isequal(a, b) && !isequal(b, c) && !isequal(a, c)
    @test hash(a) != hash(b) && hash(b) != hash(c) && hash(a) != hash(c)
end

@testset "Offset Indices" begin
    @variables k[0:3]

    @testset "i = $i" for i in 0:3
        sym = unwrap(k[i])
        @test operation(sym) === getindex
        args = arguments(sym)
        @test length(args) == 2
        @test args[1] === unwrap(k)
        @test args[2] === i
    end

    @test_throws BoundsError k[-1]
    @test_throws BoundsError k[4]
end

@testset "Arrayop sorted_arguments" begin
    @variables x[1:3] y[1:3]
    sym = unwrap(x + y)
    @test all(splat(isequal), zip(SymbolicUtils.sorted_arguments(sym), [+, x, y]))
end
