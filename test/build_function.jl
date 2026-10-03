using Symbolics, SparseArrays, LinearAlgebra, Test
using ReferenceTests
using Symbolics: value
using SymbolicUtils.Code: Assignment, DestructuredArgs, Func, LiteralExpr, NameState, Let, cse
@variables a b c1 c2 c3 d e g
oop, iip = Symbolics.build_function([sqrt(a), sin(b)], [a, b], nanmath = true)
oop = eval(oop)
@test all(isnan, @invokelatest oop([-1, Inf]))
out = [0, 0.0]
iip = eval(iip)
@invokelatest iip(out, [-1, Inf])
@test all(isnan, out)

# Multiple argument matrix
h = [a + b + c1 + c2,
     c3 + d + e + g,
     0] # uses the same number of arguments as our application
h_julia(a, b, c, d, e, g) = [a[1] + b[1] + c[1] + c[2],
                             c[3] + d[1] + e[1] + g[1],
                             0]
function h_julia!(out, a, b, c, d, e, g)
    out .= [a[1] + b[1] + c[1] + c[2], c[3] + d[1] + e[1] + g[1], 0]
end

h_str = Symbolics.build_function(h, [a], [b], [c1, c2, c3], [d], [e], [g])
h_str2 = Symbolics.build_function(h, [a], [b], [c1, c2, c3], [d], [e], [g])
@test h_str[1] == h_str2[1]
@test h_str[2] == h_str2[2]

h_oop = eval(h_str[1])
h_str_par = Symbolics.build_function(h, [a], [b], [c1, c2, c3], [d], [e], [g], parallel=Symbolics.MultithreadedForm())
h_str_3 = Symbolics.build_function(h, [a], [b], [c1, c2, c3], [d], [e], [g], iip_config = (false, true))
h_str_4 = Symbolics.build_function(h, [a], [b], [c1, c2, c3], [d], [e], [g], iip_config = (true, false))

@test contains(repr(h_str_par[1]), "schedule")
@test contains(repr(h_str_par[2]), "schedule")
h_oop_par = eval(h_str_par[1])
h_par_rgf = Symbolics.build_function(h, [a], [b], [c1, c2, c3], [d], [e], [g], parallel=Symbolics.MultithreadedForm(), expression=false)
h_ip! = eval(h_str[2])
h_ip_skip! = eval(Symbolics.build_function(h, [a], [b], [c1, c2, c3], [d], [e], [g], skipzeros=true, fillzeros=false)[2])
h_ip_skip_par! = eval(Symbolics.build_function(h, [a], [b], [c1, c2, c3], [d], [e], [g], skipzeros=true, parallel=Symbolics.MultithreadedForm(), fillzeros=false)[2])
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
@test h_oop_par(inputs...) == h_julia(inputs...)
@test h_par_rgf[1](inputs...) == h_julia(inputs...)
out_1 = similar(h, Int)
out_2 = similar(out_1)
h_ip!(out_1, inputs...)
h_julia!(out_2, inputs...)
@test_throws ArgumentError h3_oop(inputs...)
@test out_1 == out_2
h3_ip(out_1, inputs...)
@test out_1 == out_2
@test_throws ArgumentError h4_ip(out_1, inputs...)
@test h4_oop(inputs...) == h_julia(inputs...)
out_1 = similar(h, Int)
h_par_rgf[2](out_1, inputs...)
@test out_1 == out_2
fill!(out_1, 10)
h_ip_skip!(out_1, inputs...)
@test out_1[3] == 10
out_1[3] = 0
@test out_1 == out_2

fill!(out_1, 10)
h_ip_skip_par!(out_1, inputs...)
@test out_1[3] == 10
out_1[3] = 0
@test out_1 == out_2

# Multiple input matrix, some unused arguments
h_skip = [a + b + c1; c2 + c3 + g] # skip d, e
h_julia_skip(a, b, c, d, e, g) = [a[1] + b[1] + c[1]; c[2] + c[3] + g[1]]
function h_julia_skip!(out, a, b, c, d, e, g)
    out .= [a[1] + b[1] + c[1]; c[2] + c[3] + g[1]]
end

h_str_skip = Symbolics.build_function(h_skip, [a], [b], [c1, c2, c3], [], [], [g], checkbounds=true)
h_str_skip_cse = Symbolics.build_function(h_skip, [a], [b], [c1, c2, c3], [], [], [g], checkbounds=true, cse=true)
h_oop_skip = let f = eval(h_str_skip[1])
    (args...) -> @invokelatest f(args...)
end
h_ip!_skip = let f = eval(h_str_skip[2])
    (args...) -> @invokelatest f(args...)
end
h_oop_skip_cse = let f = eval(h_str_skip_cse[1])
    (args...) -> @invokelatest f(args...)
end
h_ip!_skip_cse = let f = eval(h_str_skip_cse[2])
    (args...) -> @invokelatest f(args...)
end
inputs_skip = ([1], [2], [3, 4, 5], [], [], [8])

@test h_oop_skip(inputs_skip...) == h_julia_skip(inputs_skip...) == h_oop_skip_cse(inputs_skip...)
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

# Multiple input scalar
h_scalar = a + b + c1 + c2 + c3 + d + e + g
h_julia_scalar(a, b, c, d, e, g) = a[1] + b[1] + c[1] + c[2] + c[3] + d[1] + e[1] + g[1]
h_str_scalar = Symbolics.build_function(h_scalar, [a], [b], [c1, c2, c3], [d], [e], [g])
h_str_scalar2 = Symbolics.build_function(h_scalar, [a], [b], [c1, c2, c3], [d], [e], [g])
h_str_scalar_cse = Symbolics.build_function(h_scalar, [a], [b], [c1, c2, c3], [d], [e], [g], cse=true)
@test h_str_scalar == h_str_scalar2

h_oop_scalar = let f = eval(h_str_scalar)
    (args...) -> f(args...)
end
h_oop_scalar_cse = let f = eval(h_str_scalar_cse)
    (args...) -> f(args...)
end
@test h_oop_scalar(inputs...) == h_julia_scalar(inputs...) == h_oop_scalar_cse(inputs...)

@variables z[1:100]
@variables t x(t) y(t) k
f = let _f = eval(build_function((x+y)/k, [x,y,k]))
    (args...) -> @invokelatest _f(args...)
end
@test f([1,1,2]) == 1

f = let _f = eval(build_function([(x+y)/k], [x,y,k])[1])
    (args...) -> @invokelatest _f(args...)
end
@test f([1,1,2]) == [1]

f = let _f = eval(build_function([(x+y)/k], [x,y,k])[2])
    (args...) -> @invokelatest _f(args...)
end
z = [0.0]
f(z, [1,1,2])
@test z == [1]

f = let _f = eval(build_function(sparse([1],[1], [(x+y)/k], 10,10), [x,y,k])[1])
    (args...) -> @invokelatest _f(args...)
end

@test size(f([1.,1.,2])) == (10,10)
@test f([1.,1.,2])[1,1] == 1.0
@test sum(f([1.,1.,2])) == 1.0

# Reshaped SparseMatrix optimization
let
    @variables a b c

    x = reshape(sparse([0 a 0; 0 b c]), 3, 2)
    f1,f2=build_function(x, [a,b,c], expression=Val{false})
    y = f1([1,2,3])
    @test y isa Base.ReshapedArray
    @test y.parent isa SparseMatrixCSC
    @test y.parent.rowval == x.parent.rowval
    @test y == [0 2; 0 0; 1 3]

    f1,f2=build_function(@views(x[2:3,1:2]), [a,b,c], expression=Val{false})
    y = f1([1,2,3])
    @test y isa SparseMatrixCSC
    @test y == [0 0; 1 3]
end

let # ModelingToolkit.jl#800
    @variables x
    y = sparse(1:3,1:3,x)

    f1,f2 = build_function(y,x)
    sf1, sf2 = string(f1), string(f2)
    @test !contains(sf1, "CartesianIndex")
    @test !contains(sf2, "CartesianIndex")
    @test contains(sf2, ".nzval")
end

let # Symbolics.jl#123
    ns = 6
    @variables x[1:ns]
    @variables u
    @variables M[1:36]
    @variables qd[1:6]
    output_eq = u*(qd[1]*(M[1]*qd[1] + M[7]*qd[2]))

    @test_reference "target_functions/issue123.c" build_function(output_eq, x, target=Symbolics.CTarget())
end

using Symbolics: value
using SymbolicUtils.Code: Func, toexpr
@variables t x(t)
D = Differential(t)
expr = toexpr(Func([value(D(x))], [], value(D(x))))
@test expr.args[2].args[end] == expr.args[1].args[1] # check function body and function arg
@test expr.args[2].args[end] == :(var"Differential(t, 1)(x(t))")

## Oop Arr case:
#

a = rand(4)
@variables x[1:4]
f = eval(build_function(sin.(cos.(x)), cos.(x))[1])
@test @invokelatest(f(a)) == sin.(a)

# more skipzeros
@variables x,y
f = [0, x]
f_expr = build_function(f, [x,y];skipzeros=true, expression = Val{false})

out = Vector{Float64}(undef, 2)
u = [5.0, 3.1]
@test f_expr[1](u) == [0, 5]
old = out[1]
f_expr[2](out, u)
@test out[1] === old
@test out[2] === u[1]


let # issue#136
    N = 8
    @variables x y
    A = sparse(Tridiagonal([x^i for i in 1:N-1],
                           [x^i * y^(8-i) for i in 1:N],
                           [y^i for i in 1:N-1]))

    val = Dict(x=>1, y=>2)
    B = map(A) do e
        Num(substitute(e, val))::Num
    end

    C = copy(B) - 100*I
    C_2 = copy(C);

    f = build_function(A,[x,y],parallel=Symbolics.MultithreadedForm())[2]
    g = eval(f)
    f_cse = build_function(A,[x,y],parallel=Symbolics.MultithreadedForm(),cse=true)[2]
    g_cse = eval(f_cse)

    @invokelatest g(C, [1,2])
    @test contains(repr(f), "schedule")
    @test isequal(C, B)
    @invokelatest g_cse(C_2, [1,2])
    @test isequal(C_2, B)
end


let #issue#587
    using Symbolics, SparseArrays

    N = 100 # try with N = 5 and N = 100
    _S = sprand(N, N, 0.1)
    _Q = Array(sprand(N, N, 0.1))

    F(z) = [
            collect(_S * z)
            collect(_Q * z.^2)
           ]

    Symbolics.@variables z[1:N]

    sj = Symbolics.sparsejacobian(F(z), z)

    f_expr = build_function(sj, z)
    myf = eval(first(f_expr))
    J = @invokelatest myf(rand(N))

    @test typeof(J) <: SparseMatrixCSC
    @test nnz(J) == nnz(sj)
end

let # Symbolics.jl#2006
    @variables a b
    # empty sparse output: every entry is structurally zero
    H_empty = Symbolics.sparsehessian(a + b, [a, b])
    @test isempty(H_empty.nzval)
    f_empty = eval(build_function(H_empty, [a, b])[1])
    out_empty = @invokelatest f_empty([1.0, 2.0])
    # must support arithmetic against a numeric matrix
    @test out_empty ≈ zeros(2, 2)
    # eltype must be concrete, not `Any`
    @test eltype(out_empty) <: Number
    # same concrete eltype as the non-empty sparse output for the same inputs
    H_const = Symbolics.sparsehessian(a^2 + b, [a, b])
    f_const = eval(build_function(H_const, [a, b])[1])
    out_const = @invokelatest f_const([1.0, 2.0])
    @test eltype(out_empty) == eltype(out_const)
end

# test header wrapping of scalar build function
let
    @variables x p t
    ex = t + p * x^2
    integrator = gensym(:MTKIntegrator)
    header = expr -> let integrator = integrator
        Func([expr.args[1], expr.args[2], DestructuredArgs(expr.args[3:end], integrator,
                                                           inds = [:p])], [], expr.body)
    end
    f = build_function(ex, [value(x)], value(t), [value(p)]; expression = Val{false},
                                                                wrap_code = header)
    p = (a = 10, p = [2])
    @test f([3], 1, p) == 19
end

let #658
    using Symbolics
    @variables a, X1[1:3], X2[1:3]
    k = eval(build_function(a * X1 + X2, X1, X2, a)[1])
    @test @invokelatest(k(ones(3), ones(3), 1.5)) == [2.5, 2.5, 2.5]
end

@testset "ArrayOp codegen" begin
    @variables x[1:2]
    T = value(x .^ 2)
    @test_nowarn toexpr(T, NameState())
end

@testset "`similarto` keyword argument" begin
    @variables x[1:2]
    T = collect(value(x .^ 2))
    fn = build_function(T, collect(x); expression = false)[1]
    @test_throws MethodError fn((1.0, 2.0))
    fn = build_function(T, collect(x); similarto = Array, expression = false)[1]
    @test fn((1.0, 2.0)) ≈ [1.0, 4.0]
end

@testset "`build_function` with array symbolics" begin
    @variables x[1:4]
    for var in [x[1:2], x[1:2] .+ 0.0, Symbolics.unwrap(x[1:2])]
        foop, fiip = build_function(var[1:2], x; expression = false)
        @test foop(ones(4)) ≈ ones(2)
        buf = zeros(2)
        fiip(buf, ones(4))
        @test buf ≈ ones(2)
    end
end

@testset "cse with arrayops" begin
    @variables x[1:3] y f(..)
    t = x .+ y
    t = t .* f(t)
    res = cse(value(t))
    @test res isa Let
    @test !isempty(res.pairs)
end

@testset "`CallWithMetadata` in `DestructuredArgs` with `create_bindings = false`" begin
    @variables x f(..)
    fn = build_function(f(x), DestructuredArgs([f]; create_bindings = false), x; expression = Val{false})
    @test fn([isodd], 3)
end

@testset "iip_config with RGF" begin
    @variables a b
    oop, iip = build_function([a + b, a - b], a, b; iip_config = (false, false), expression = Val{false})
    @test_throws ArgumentError oop(1, 2)
    @test_throws ArgumentError iip(ones(2), 1, 2)

    @variables a[1:2]
    oop, iip = build_function(a .* 2, a; iip_config = (false, false), expression = Val{false})
    @test_throws ArgumentError oop(ones(2))
    @test_throws ArgumentError iip(ones(2), ones(2))
end

@testset "unwrapping/CSE in array of symbolics codegen" begin
    @variables a b
    oop, _ = build_function([a^2 + b^2, a^2 + b^2], a, b; expression = Val{true}, cse = true)

    function find_create_array(expr)
        while expr isa Expr && (!Meta.isexpr(expr, :call) || expr.args[1] != SymbolicUtils.Code.create_array)
            expr = expr.args[end]
        end
        return expr
    end

    expr = find_create_array(oop)
    # CSE works, we just need to test that it's happening and OOP is the easiest way to do it
    @test Meta.isexpr(expr, :call) && expr.args[1] == SymbolicUtils.Code.create_array &&
          expr.args[end] isa Symbol && expr.args[end-1] isa Symbol
end

@testset "CSE with operators" begin
    @variables t x(t)
    D = Differential(t)
    f = build_function(x + D(x), [x, D(x)]; cse = true, expression = Val{false})
    @test f([1, 2]) == 3
end

@testset "`build_function` with `UpperTriangular`" begin
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

    out, fjac_upper_expr = build_function(up_J, u; skipzeros = true, expression = false)
    Jtmp = UpperTriangular(zeros(2, 2))
    utmp = rand(2)
    @test_nowarn fjac_upper_expr(Jtmp, utmp)
    @test Jtmp[3] == utmp[2]
end

@testset "MultithreadedForm expressions can be written to a file and included" begin
    @variables x y
    A = [
        x^2 + y 0 2x
        0 0 2y
        y^2 + x 0 0
    ]
    u = [1.0, 2.0]
    expected = @invokelatest eval(build_function(A, [x, y])[1])(u)
    oop_ex, iip_ex = build_function(A, [x, y]; parallel = Symbolics.MultithreadedForm())
    mktempdir() do dir
        oop_path = joinpath(dir, "f_oop.jl")
        iip_path = joinpath(dir, "f_iip.jl")
        write(oop_path, string(oop_ex))
        write(iip_path, string(iip_ex))
        f_oop = include(oop_path)
        f_iip = include(iip_path)
        @test @invokelatest(f_oop(u)) == expected
        out = zeros(3, 3)
        @invokelatest f_iip(out, u)
        @test out == expected
    end
end

# `-tN,0` keeps the main thread in the default pool on 1.12+; older Julia rejects `,0`.
const MT_TEST_THREADS = VERSION >= v"1.12" ? "4,0" : "4"

@testset "MultithreadedForm in-place function from Threads.@spawn" begin
    # A race between the generated tasks can crash the whole Julia process, so
    # the calls have to run in a subprocess with real worker threads.
    script = joinpath(mktempdir(), "mt_iip_spawn.jl")
    write(
        script, """
        using Symbolics
        @variables x y
        N = 8
        A = Num[x^i + y^j for i in 1:N, j in 1:N]
        u = [1.0, 2.0]
        _, f_serial = build_function(A, [x, y]; parallel = Symbolics.SerialForm(), expression = Val(false))
        _, f_par = build_function(A, [x, y]; parallel = Symbolics.MultithreadedForm(2, 4), expression = Val(false))
        ref = zeros(N, N)
        f_serial(ref, u)
        f_par(zeros(N, N), u)
        outs = [zeros(N, N) for _ in 1:200]
        ok = all(fetch, map(1:200) do i
            Threads.@spawn begin
                f_par(outs[i], u)
                outs[i] == ref
            end
        end)
        println(ok ? "ALL_CORRECT" : "MISMATCH")
        exit(ok ? 0 : 1)
        """
    )
    cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) -t$(MT_TEST_THREADS) $script`
    @test success(pipeline(cmd; stdout = stdout, stderr = stderr))
end

@testset "MultithreadedForm RuntimeGeneratedFunction with three or more arguments" begin
    # Julia 1.10/1.11 segfault when this generated code calls opaque closures, so the
    # calls run in a subprocess.
    script = joinpath(mktempdir(), "mt_rgf_nargs.jl")
    write(
        script, """
        using Symbolics
        @variables a b c d
        h = [a + b + c, c + d, a * d, 0]
        args = ([a], [b], [c], [d])
        inputs = ([1], [2], [3], [4])
        expected = [6, 7, 4, 0]
        for nt in (Symbolics.MultithreadedForm(), Symbolics.MultithreadedForm(2, 4))
            f_oop, f_iip = build_function(h, args...; parallel = nt, expression = Val(false))
            f_oop(inputs...) == expected || exit(1)
            out = zeros(Int, 4)
            f_iip(out, inputs...)
            out == expected || exit(1)
        end
        println("ALL_CORRECT")
        """
    )
    for threads in (1, MT_TEST_THREADS)
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) -t$(threads) $script`
        @test success(pipeline(cmd; stdout = stdout, stderr = stderr))
    end
end

@testset "Shard postprocessing scope" begin
    @variables x y z
    for parallel in (Symbolics.SerialForm(), Symbolics.ShardedForm(1, 2), Symbolics.MultithreadedForm(1, 2)),
            expression in (Val{true}, Val{false}), cse in (false, true),
            sequential in (false, true), sparse_output in (false, true)
        ex, expected = sparse_output ?
            (sparse([z + 1 0; z + 2 z]), sparse([7.0 0; 8.0 6.0])) :
            ([z + 1, z + 2, x, y], [7.0, 8.0, 2.0, 3.0])
        wrap = if sequential
            b -> Let([Assignment(z, x * y), Assignment(:result, b)], :result, false)
        else
            b -> Let([Assignment(z, x * y)], b, false)
        end
        @testset "$parallel $expression cse=$cse sequential=$sequential sparse=$sparse_output" begin
            fs = build_function(ex, [x, y]; parallel, expression, cse, postprocess_fbody = wrap)
            f, g = expression == Val{true} ? eval.(fs) : fs
            @test Base.invokelatest(f, [2.0, 3.0]) == expected
            out = copy(expected)
            fill!(out, 0)
            Base.invokelatest(g, out, [2.0, 3.0])
            @test out == expected
        end
    end
end

@testset "Whole-result shard postprocessing" begin
    @variables x z
    for parallel in (Symbolics.SerialForm(), Symbolics.ShardedForm(1, 2), Symbolics.MultithreadedForm(1, 2)),
            expression in (Val{true}, Val{false})
        fs = build_function(
            [x, 2x, 3x, 4x], x;
            parallel, expression, iip_config = (true, false),
            postprocess_fbody = b -> LiteralExpr(:(reverse($b)))
        )
        f = expression == Val{true} ? eval(first(fs)) : first(fs)
        @test Base.invokelatest(f, 1) == [4, 3, 2, 1]

        @testset "ordinary closure $parallel $expression" begin
            wrap = b -> Let([Assignment(z, 3)], LiteralExpr(:($(Func([], [], b))())), false)
            fs = build_function([x, 2x, 3x, 4x], x; parallel, expression, postprocess_fbody = wrap)
            f, g = expression == Val{true} ? eval.(fs) : fs
            @test Base.invokelatest(f, 2) == [2, 4, 6, 8]
            out = zeros(Int, 4)
            Base.invokelatest(g, out, 2)
            @test out == [2, 4, 6, 8]
        end

        @testset "indexed mutation $parallel $expression" begin
            counter = Ref(0)
            wrap = b -> Let(
                [
                    Assignment(LiteralExpr(:($counter[])), LiteralExpr(:($counter[] + 1))),
                ], b, false
            )
            fs = build_function([x, 2x, 3x, 4x], x; parallel, expression, postprocess_fbody = wrap)
            f, g = expression == Val{true} ? eval.(fs) : fs
            @test Base.invokelatest(f, 2) == [2, 4, 6, 8]
            @test counter[] == 1
            out = zeros(Int, 4)
            Base.invokelatest(g, out, 2)
            @test out == [2, 4, 6, 8]
            @test counter[] == 2
        end
    end
end

@testset "Shard bindings and input destructuring" begin
    @variables x y z
    for parallel in (Symbolics.SerialForm(), Symbolics.ShardedForm(1, 2), Symbolics.MultithreadedForm(1, 2))
        @testset "input mutation $parallel" begin
            wrap = b -> Let([Assignment(x, 3x)], b, false)
            f, g = build_function(
                [x, y, 2x, 2y], [x, y];
                parallel, expression = Val{false}, postprocess_fbody = wrap
            )
            @test f([2, 3]) == [6, 3, 12, 6]
            out = zeros(Int, 4)
            g(out, [2, 3])
            @test out == [6, 3, 12, 6]
        end
        for expression in (Val{false}, Val(false), false), wrap in (
                    b -> Let([Assignment(z, x * y)], b, true),
                    b -> Let([Assignment(z, 1)], Let([Assignment(z, x * y)], b, true), true),
                    b -> Let([Assignment(DestructuredArgs([z], :values), LiteralExpr(:([6])))], b, false),
                    b -> Let([DestructuredArgs([z], LiteralExpr(:([6])))], b, false),
                    b -> Let([Assignment(LiteralExpr(:((z,))), LiteralExpr(:((6,))))], b, false),
                    b -> Let([Assignment(DestructuredArgs((z,), :values), LiteralExpr(:([6])))], b, false),
                    b -> Let([Assignment(LiteralExpr(:((_, z))), LiteralExpr(:((1, 6))))], b, false),
                )
            f, g = build_function(
                [z + 1, z + 2, x, y], [x, y];
                parallel, expression, postprocess_fbody = wrap
            )
            @test f([2, 3]) == [7, 8, 2, 3]
            out = zeros(Int, 4)
            g(out, [2, 3])
            @test out == [7, 8, 2, 3]
        end
    end
end
