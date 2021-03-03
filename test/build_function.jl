using Symbolics, SparseArrays, Test
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
h_oop = eval(h_str[1])
h_ip! = eval(h_str[2])
h_ip_skip! = eval(Symbolics.build_function(h, [a], [b], [c1, c2, c3], [d], [e], [g], skipzeros=true, fillzeros=false)[2])
inputs = ([1], [2], [3, 4, 5], [6], [7], [8])

@test h_oop(inputs...) == h_julia(inputs...)
out_1 = similar(h, Int)
out_2 = similar(out_1)
h_ip!(out_1, inputs...)
h_julia!(out_2, inputs...)
@test out_1 == out_2
fill!(out_1, 10)
h_ip_skip!(out_1, inputs...)
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
h_oop_scalar = eval(h_str_scalar)
@test h_oop_scalar(inputs...) == h_julia_scalar(inputs...)

@variables z[1:100]
@test isequal(simplify(Symbolics.unflatten_long_ops(sum(z))),
              simplify(sum(z)))

@test isequal(simplify(Symbolics.unflatten_long_ops(prod(z))),
              simplify(prod(z)))

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

# Matrix CTarget test
let 
    @variables x[1:4] y[1:4] z[1:4]
    expression = hcat(x,y,z) # Results in [x1 y1 z1; x2 y2 z2; ...]
    variables  = vcat(x,y,z) # Results in [x1, x2, x3, x4, y1, ...]
    cfunc      = build_function(expression, variables; target = Symbolics.CTarget(), expression = Val{true})

    # Generated function should be out[0] = in[0], out[1] = in[1], out[2] = in[2], etc.
    @test cfunc == "void diffeqf(double* du, double* RHS1) { \n  du[0] = RHS1[0];\n  du[1] = RHS1[1];\n  du[2] = RHS1[2];\n  du[3] = RHS1[3];\n  du[4] = RHS1[4];\n  du[5] = RHS1[5];\n  du[6] = RHS1[6];\n  du[7] = RHS1[7];\n  du[8] = RHS1[8];\n  du[9] = RHS1[9];\n  du[10] = RHS1[10];\n  du[11] = RHS1[11]; \n}\n"
end

# Scalar CTarget test
let 
    @variables x y z
    expression = x + y + z
    cfunc = build_function(expression, [x], [y], [z]; target = Symbolics.CTarget(), expression = Val{true})

    # Generated function should be out[0] = in1[0] + in2[0] + in3[0]
    @test cfunc == "void diffeqf(double* du, double* RHS1, double* RHS2, double* RHS3) { \n  du[0] = RHS1[0] + RHS2[0] + RHS3[0]; \n}\n"
end

# Matrix StanTarget test
let 
    @variables x[1:4] y[1:4] z[1:4]
    expression = hcat(x,y,z) # Results in [x1 y1 z1; x2 y2 z2; ...]
    variables  = vcat(x,y,z) # Results in [x1, x2, x3, x4, y1, ...]
    sfunc      = build_function(expression, vcat(x,y), z, []; target = Symbolics.StanTarget(), expression = Val{true})

    # Generated function should be out[0] = v[0], out[1] = v[1], out[2] = v[2], ..., out[8] = p[8], etc.
    @test sfunc == "real[] diffeqf(real Any[],real[] internal_var___u,real[] internal_var___p,real[] x_r,int[] x_i) {\n  real internal_var___du[12];\n\n  internal_var___du[0] = internal_var___u[0];\n  internal_var___du[1] = internal_var___u[1];\n  internal_var___du[2] = internal_var___u[2];\n  internal_var___du[3] = internal_var___u[3];\n  internal_var___du[4] = internal_var___u[4];\n  internal_var___du[5] = internal_var___u[5];\n  internal_var___du[6] = internal_var___u[6];\n  internal_var___du[7] = internal_var___u[7];\n  internal_var___du[8] = internal_var___p[0];\n  internal_var___du[9] = internal_var___p[1];\n  internal_var___du[10] = internal_var___p[2];\n  internal_var___du[11] = internal_var___p[3]; \n\n  return internal_var___du;\n}\n"
end

# Scalar StanTarget test
let 
    @variables t x(t) y(t) z(t)
    expression = x + y + z
    sfunc = build_function(expression, vcat(x,y), [z], t; target = Symbolics.StanTarget(), expression = Val{true})

    # Generated function should be out[0] = v[0] + p[0] + t[0]
    @test sfunc == "real[] diffeqf(real t,real[] internal_var___u,real[] internal_var___p,real[] x_r,int[] x_i) {\n  real internal_var___du[1];\n\n  internal_var___du[0] = internal_var___u[0] + internal_var___u[1] + internal_var___p[0]; \n\n  return internal_var___du;\n}\n"
end

# Matrix MATLABTarget test
let 
    @variables x[1:4] y[1:4] z[1:4]
    expression = hcat(x,y,z) # Results in [x1 y1 z1; x2 y2 z2; ...]
    variables  = vcat(x,y,z) # Results in [x1, x2, x3, x4, y1, ...]
    mfunc      = build_function(expression, vcat(x,y,z); target = Symbolics.MATLABTarget(), expression = Val{true})

    # Generated function should be of the same form as expression
    @test mfunc == "diffeqf = @(t,internal_var___u) [\n  internal_var___u(1), internal_var___u(5), internal_var___u(9);\n  internal_var___u(2), internal_var___u(6), internal_var___u(10);\n  internal_var___u(3), internal_var___u(7), internal_var___u(11);\n  internal_var___u(4), internal_var___u(8), internal_var___u(12);\n];"
end

# Scalar MATLABTarget test
let 
    @variables x y z
    expression = x + y + z
    mfunc = build_function(expression, vcat(x,y,z); target = Symbolics.MATLABTarget(), expression = Val{true}, allowscalar = true)

    # Generated function should be out[0] = v[0] + p[0] + t[0]
    @test mfunc == "diffeqf = @(t,internal_var___u)  internal_var___u(1) + internal_var___u(2) + internal_var___u(3);\n;"
end

