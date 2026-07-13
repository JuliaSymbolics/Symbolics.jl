using Symbolics
using SymbolicUtils
using SymbolicUtils.Code
using SymbolicUtils: unwrap
using SparseArrays
using LinearAlgebra
using RuntimeGeneratedFunctions

RuntimeGeneratedFunctions.init(@__MODULE__)

@variables a b c1 c2 c3 d e g

ir = IRStructure{SymReal}()

@testset "NaNMath" begin
    oop, iip = Symbolics.codegen_function(ir, [sqrt(a), sin(b)], [[a, b]]; nanmath = true, sort_addmul = true)
    oop = eval(oop)
    @test all(isnan, @invokelatest oop([-1, Inf]))
    out = [0, 0.0]
    iip = eval(iip)
    @invokelatest iip(out, [-1, Inf])
    @test all(isnan, out)
end

@testset "Multiple arguments" begin
    h = [
        a + b + c1 + c2,
        c3 + d + e + g,
        0,
    ] # uses the same number of arguments as our application
    h_julia(a, b, c, d, e, g) = [
        a[1] + b[1] + c[1] + c[2],
        c[3] + d[1] + e[1] + g[1],
        0,
    ]
    function h_julia!(out, a, b, c, d, e, g)
        out .= [a[1] + b[1] + c[1] + c[2], c[3] + d[1] + e[1] + g[1], 0]
    end

    h_str = Symbolics.codegen_function(ir, h, [[a], [b], [c1, c2, c3], [d], [e], [g]]; sort_addmul = true)
    h_str2 = Symbolics.codegen_function(ir, h, [[a], [b], [c1, c2, c3], [d], [e], [g]]; sort_addmul = true)
    @test h_str[1] == h_str2[1]
    @test h_str[2] == h_str2[2]

    h_oop = eval(h_str[1])
    h_str_3 = Symbolics.codegen_function(ir, h, [[a], [b], [c1, c2, c3], [d], [e], [g]]; iip_config = (false, true), sort_addmul = true)
    h_str_4 = Symbolics.codegen_function(ir, h, [[a], [b], [c1, c2, c3], [d], [e], [g]]; iip_config = (true, false), sort_addmul = true)

    h_ip! = eval(h_str[2])
    h_ip_skip! = eval(Symbolics.codegen_function(ir, h, [[a], [b], [c1, c2, c3], [d], [e], [g]], skipzeros = true, sort_addmul = true)[2])

    h3_oop = let f = eval(h_str_3[1])
        (args...) -> @invokelatest f(args...)
    end
    h3_ip = let f = eval(h_str_3[2])
        (args...) -> @invokelatest f(args...)
    end
    h4_oop = let f = eval(h_str_4[1])
        (args...) -> @invokelatest f(args...)
    end
    h4_ip = let f = eval(h_str_4[2])
        (args...) -> @invokelatest f(args...)
    end
    inputs = ([1], [2], [3, 4, 5], [6], [7], [8])

    @test h_oop(inputs...) == h_julia(inputs...)
    out_1 = similar(h, Int)
    out_2 = similar(out_1)
    h_ip!(out_1, inputs...)
    h_julia!(out_2, inputs...)
    @test_throws Symbolics.FunctionUnimplementedError h3_oop(inputs...)
    @test out_1 == out_2
    h3_ip(out_1, inputs...)
    @test out_1 == out_2
    @test_throws Symbolics.FunctionUnimplementedError h4_ip(out_1, inputs...)
    @test h4_oop(inputs...) == h_julia(inputs...)
    out_1 = similar(h, Int)
    fill!(out_1, 10)
    h_ip_skip!(out_1, inputs...)
    @test out_1[3] == 10
    out_1[3] = 0
    @test out_1 == out_2
end

@testset "With unused arguments" begin
    h_skip = [a + b + c1; c2 + c3 + g] # skip d, e
    h_julia_skip(a, b, c, d, e, g) = [a[1] + b[1] + c[1]; c[2] + c[3] + g[1]]
    function h_julia_skip!(out, a, b, c, d, e, g)
        out .= [a[1] + b[1] + c[1]; c[2] + c[3] + g[1]]
    end

    h_str_skip = Symbolics.codegen_function(ir, h_skip, [[a], [b], [c1, c2, c3], [], [], [g]]; checkbounds = true, sort_addmul = true)
    h_oop_skip = let f = eval(h_str_skip[1])
        (args...) -> @invokelatest f(args...)
    end
    h_ip!_skip = let f = eval(h_str_skip[2])
        (args...) -> @invokelatest f(args...)
    end
    inputs_skip = ([1], [2], [3, 4, 5], [], [], [8])

    @test h_oop_skip(inputs_skip...) == h_julia_skip(inputs_skip...)
    out_1_skip = Array{Int64}(undef, 2)
    out_2_skip = similar(out_1_skip)
    h_ip!_skip(out_1_skip, inputs_skip...)
    h_julia_skip!(out_2_skip, inputs_skip...)
    @test out_1_skip == out_2_skip

    # Same as above, except test ability to call with non-matrix arguments (i.e., for `nt`)
    inputs_skip_2 = ([1], [2], [3, 4, 5], [], (a = 1, b = 2), [8])
    @test h_oop_skip(inputs_skip_2...) == h_julia_skip(inputs_skip_2...)
    out_1_skip_2 = Array{Int64}(undef, 2)
    out_2_skip_2 = similar(out_1_skip_2)
    h_ip!_skip(out_1_skip_2, inputs_skip_2...)
    h_julia_skip!(out_2_skip_2, inputs_skip_2...)
    @test out_1_skip_2 == out_2_skip_2
end

@testset "Multiple input scalar result" begin
    h_scalar = a + b + c1 + c2 + c3 + d + e + g
    h_julia_scalar(a, b, c, d, e, g) = a[1] + b[1] + c[1] + c[2] + c[3] + d[1] + e[1] + g[1]
    h_str_scalar, _ = Symbolics.codegen_function(ir, h_scalar, [[a], [b], [c1, c2, c3], [d], [e], [g]]; sort_addmul = true)
    h_str_scalar2, _ = Symbolics.codegen_function(ir, h_scalar, [[a], [b], [c1, c2, c3], [d], [e], [g]]; sort_addmul = true)
    @test h_str_scalar == h_str_scalar2

    h_oop_scalar = let f = eval(h_str_scalar)
        (args...) -> f(args...)
    end
    inputs = ([1], [2], [3, 4, 5], [6], [7], [8])
    @test h_oop_scalar(inputs...) == h_julia_scalar(inputs...)
end

@testset "Dependent variable arguments" begin
    @variables t x(t) y(t) k
    f = let _f = eval(Symbolics.codegen_function(ir, (x + y) / k, [[x, y, k]]; sort_addmul = true)[1])
        (args...) -> @invokelatest _f(args...)
    end
    @test f([1, 1, 2]) == 1

    f = let _f = eval(Symbolics.codegen_function(ir, [(x + y) / k], [[x, y, k]]; sort_addmul = true)[1])
        (args...) -> @invokelatest _f(args...)
    end
    @test f([1, 1, 2]) == [1]

    f = let _f = eval(Symbolics.codegen_function(ir, [(x + y) / k], [[x, y, k]]; sort_addmul = true)[2])
        (args...) -> @invokelatest _f(args...)
    end
    z = [0.0]
    f(z, [1, 1, 2])
    @test z == [1]

    f = let _f = eval(Symbolics.codegen_function(ir, sparse([1], [1], [(x + y) / k], 10, 10), [[x, y, k]]; sort_addmul = true)[1])
        (args...) -> @invokelatest _f(args...)
    end

    @test size(f([1.0, 1.0, 2])) == (10, 10)
    @test f([1.0, 1.0, 2])[1, 1] == 1.0
    @test sum(f([1.0, 1.0, 2])) == 1.0
end

@testset "Reshaped SparseMatrix" begin
    @variables a b c

    x = reshape(sparse([0 a 0; 0 b c]), 3, 2)
    f1, f2 = Symbolics.codegen_function(ir, x, [[a, b, c]]; sort_addmul = true)
    f1 = @RuntimeGeneratedFunction(f1)
    f2 = @RuntimeGeneratedFunction(f2)
    y = f1([1, 2, 3])
    @test y isa Base.ReshapedArray
    @test y.parent isa SparseMatrixCSC
    @test y.parent.rowval == x.parent.rowval
    @test y == [0 2; 0 0; 1 3]

    f1, f2 = Symbolics.codegen_function(ir, @views(x[2:3, 1:2]), [[a, b, c]]; sort_addmul = true)
    f1 = @RuntimeGeneratedFunction(f1)
    f2 = @RuntimeGeneratedFunction(f2)
    y = f1([1, 2, 3])
    @test y isa SparseMatrixCSC
    @test y == [0 0; 1 3]
end

@testset "ModelingToolkit.jl#800" begin
    @variables x
    y = sparse(1:3, 1:3, x)

    f1, f2 = Symbolics.codegen_function(ir, y, [x]; sort_addmul = true)
    sf1, sf2 = string(f1), string(f2)
    @test !contains(sf1, "CartesianIndex")
    @test !contains(sf2, "CartesianIndex")
    @test contains(sf2, ".nzval")
end

@testset "Issue#587" begin
    N = 100 # try with N = 5 and N = 100
    _S = sprand(N, N, 0.1)
    _Q = Array(sprand(N, N, 0.1))

    F(z) = [
        collect(_S * z)
        collect(_Q * z .^ 2)
    ]

    Symbolics.@variables z[1:N]

    sj = Symbolics.sparsejacobian(F(z), z)

    f_expr = Symbolics.codegen_function(ir, sj, [z]; sort_addmul = true)
    myf = eval(first(f_expr))
    J = @invokelatest myf(rand(N))

    @test typeof(J) <: SparseMatrixCSC
    @test nnz(J) == nnz(sj)
end

@testset "`wrap_code`" begin
    @variables x p t
    ex = t + p * x^2
    integrator = gensym(:MTKIntegrator)
    header = expr -> let integrator = integrator
        Func(
            [
                expr.args[1], expr.args[2], DestructuredArgs(
                    expr.args[3:end], integrator,
                    inds = [:p]
                ),
            ], [], expr.body
        )
    end
    f, _ = Symbolics.codegen_function(
        ir, ex, [[unwrap(x)], unwrap(t), [unwrap(p)]];
        wrap_code = (header, identity), sort_addmul = true
    )
    f = @RuntimeGeneratedFunction(f)
    p = (a = 10, p = [2])
    @test f([3], 1, p) == 19
end

@testset "Issue#658" begin
    @variables a, X1[1:3], X2[1:3]
    k = eval(Symbolics.codegen_function(ir, a * X1 + X2, [X1, X2, a]; sort_addmul = true)[1])
    @test @invokelatest(k(ones(3), ones(3), 1.5)) == [2.5, 2.5, 2.5]
end

@testset "`similarto`" begin
    @variables x[1:2]
    T = collect(unwrap(x .^ 2))
    fn = @RuntimeGeneratedFunction(Symbolics.codegen_function(ir, T, [collect(x)]; sort_addmul = true)[1])
    @test_throws MethodError fn((1.0, 2.0))
    fn = @RuntimeGeneratedFunction(Symbolics.codegen_function(ir, T, [collect(x)]; similarto = Array, sort_addmul = true)[1])
    @test fn((1.0, 2.0)) ≈ [1.0, 4.0]
end

@testset "`codegen_function` with `UpperTriangular`" begin
    function f_test(J,u)
        J[1,1] = u[1]
        J[1,2] = u[2]
        J[2,1] = -u[1]
        J[2,2] = -u[2]
        return nothing
    end

    @variables u[1:2]
    J = fill!(Array{Num}(undef, 2, 2), 0)
    f_test(J, u)
    up_J = UpperTriangular(J - Diagonal(J))

    out, fjac_upper_expr = Symbolics.codegen_function(ir, up_J, [u]; skipzeros = true, sort_addmul = true)
    fjac_upper_expr = @RuntimeGeneratedFunction(fjac_upper_expr)
    Jtmp = UpperTriangular(zeros(2, 2))
    utmp = rand(2)
    @test_nowarn fjac_upper_expr(Jtmp, utmp)
    @test Jtmp[3] == utmp[2]
end

@testset "Repeated codegen produces identical expressions" begin
    @variables t x(t) fn(..) y(t)
    expr = [x + fn(y + 2t), fn(x + 3sin(t)) * y]
    args = [[x, y], [fn], t]
    oopexprs = Expr[]
    iipexprs = Expr[]
    for i in 1:10
        foop, fiip = Symbolics.codegen_function(ir, expr, args; sort_addmul = true)
        push!(oopexprs, foop)
        push!(iipexprs, fiip)
    end
    @test allequal(oopexprs)
    @test allequal(iipexprs)
end


@testset "`CodegenFunctionOptions` struct interface" begin
    @variables a b c1 c2 c3 d e g
    @test (@doc Symbolics.CodegenFunctionOptions) !== nothing
    @test (@doc Symbolics.codegen_function) !== nothing
    @static if VERSION >= v"1.11"
        @test Base.ispublic(Symbolics, :CodegenFunctionOptions)
        @test Base.ispublic(Symbolics, :codegen_function)
    end
    # `CodegenFunctionOptions` is a single concrete type regardless of which options are set, so
    # functions that thread it through do not recompile per option-combination.
    @test isconcretetype(Symbolics.CodegenFunctionOptions)
    o1 = Symbolics.CodegenFunctionOptions(; nanmath = true)
    o2 = Symbolics.CodegenFunctionOptions(;
        nanmath = false, wrap_code = (x -> x, x -> x), checkbounds = true,
        iip_config = (true, false), sort_addmul = true, skipzeros = true, optimize = nothing)
    @test typeof(o1) === typeof(o2) === Symbolics.CodegenFunctionOptions

    # scalar/vector path: keyword form and struct form are identical
    kw = Symbolics.codegen_function(ir, [sqrt(a), sin(b)], [[a, b]]; nanmath = true, sort_addmul = true)
    st = Symbolics.codegen_function(ir, [sqrt(a), sin(b)], [[a, b]],
        Symbolics.CodegenFunctionOptions(; nanmath = true, sort_addmul = true))
    @test string(kw) == string(st)

    # array path: keyword form and struct form are identical
    h = [a + b, c1 * c2, d / e]
    args = [[a], [b], [c1, c2, c3], [d], [e], [g]]
    kw2 = Symbolics.codegen_function(ir, h, args; skipzeros = true, checkbounds = true, sort_addmul = true)
    st2 = Symbolics.codegen_function(ir, h, args,
        Symbolics.CodegenFunctionOptions(; skipzeros = true, checkbounds = true, sort_addmul = true))
    @test string(kw2) == string(st2)

    # unknown keyword arguments are dropped (backwards compatible with the old `kwargs...` sink)
    with_bogus = Symbolics.codegen_function(ir, [sqrt(a)], [[a]]; nanmath = true, sort_addmul = true, bogus = 42)
    without = Symbolics.codegen_function(ir, [sqrt(a)], [[a]]; nanmath = true, sort_addmul = true)
    @test string(with_bogus) == string(without)
end
