using Symbolics, SparseArrays, LinearAlgebra, Test
using ReferenceTests
@variables a b c1 c2 c3 d e g

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

@test contains(repr(h_str_par[1]), "schedule")
@test contains(repr(h_str_par[2]), "schedule")
h_oop_par = eval(h_str_par[1])
h_par_rgf = Symbolics.build_function(h, [a], [b], [c1, c2, c3], [d], [e], [g], parallel=Symbolics.MultithreadedForm(), expression=false)
h_ip! = eval(h_str[2])
h_ip_skip! = eval(Symbolics.build_function(h, [a], [b], [c1, c2, c3], [d], [e], [g], skipzeros=true, fillzeros=false)[2])
h_ip_skip_par! = eval(Symbolics.build_function(h, [a], [b], [c1, c2, c3], [d], [e], [g], skipzeros=true, parallel=Symbolics.MultithreadedForm(), fillzeros=false)[2])
inputs = ([1], [2], [3, 4, 5], [6], [7], [8])

@test h_oop(inputs...) == h_julia(inputs...)
@test h_oop_par(inputs...) == h_julia(inputs...)
@test h_par_rgf[1](inputs...) == h_julia(inputs...)
out_1 = similar(h, Int)
out_2 = similar(out_1)
h_ip!(out_1, inputs...)
h_julia!(out_2, inputs...)
@test out_1 == out_2
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
h_oop_skip = eval(h_str_skip[1])
h_ip!_skip = eval(h_str_skip[2])
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

# Multiple input scalar
h_scalar = a + b + c1 + c2 + c3 + d + e + g
h_julia_scalar(a, b, c, d, e, g) = a[1] + b[1] + c[1] + c[2] + c[3] + d[1] + e[1] + g[1]
h_str_scalar = Symbolics.build_function(h_scalar, [a], [b], [c1, c2, c3], [d], [e], [g])
h_str_scalar2 = Symbolics.build_function(h_scalar, [a], [b], [c1, c2, c3], [d], [e], [g])
@test h_str_scalar == h_str_scalar2

h_oop_scalar = eval(h_str_scalar)
@test h_oop_scalar(inputs...) == h_julia_scalar(inputs...)

@variables z[1:100]
@variables t x(t) y(t) k
f = eval(build_function((x+y)/k, [x,y,k]))
@test f([1,1,2]) == 1

f = eval(build_function([(x+y)/k], [x,y,k])[1])
@test f([1,1,2]) == [1]

f = eval(build_function([(x+y)/k], [x,y,k])[2])
z = [0.0]
f(z, [1,1,2])
@test z == [1]

f = eval(build_function(sparse([1],[1], [(x+y)/k], 10,10), [x,y,k])[1])

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
@test expr.args[2].args[end] == :(var"Differential(t)(x(t))")

## Oop Arr case:
#

a = rand(4)
@variables x[1:4]
@test eval(build_function(sin.(cos.(x)), cos.(x)))(a) == sin.(a)

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
        Num(substitute(e, val))
    end

    C = copy(B) - 100*I


    f = build_function(A,[x,y],parallel=Symbolics.MultithreadedForm())[2]
    g = eval(f)

    g(C, [1,2])
    @test contains(repr(f), "schedule")
    @test isequal(C, B)
end
