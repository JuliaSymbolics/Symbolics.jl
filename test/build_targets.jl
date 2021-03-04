using Symbolics, Test

@variables t a x(t) y(t)
expr = [a*x - x*y,-3y + x*y]
@test Symbolics.build_function(expr,[x,y],[a],t,target = Symbolics.StanTarget()) ==
    """
    real[] diffeqf(real t,real[] internal_var___u,real[] internal_var___p,real[] x_r,int[] x_i) {
      real internal_var___du[2];
      internal_var___du[1] = internal_var___p[1] * internal_var___u[1] + -1 * internal_var___u[1] * internal_var___u[2];
      internal_var___du[2] = internal_var___u[1] * internal_var___u[2] + -3 * internal_var___u[2];
      return internal_var___du;
    }
    """

@test Symbolics.build_function(expr,[x,y],[a],t,target = Symbolics.CTarget(),
                                     lhsname=:internal_var___du,
                                     rhsnames=[:internal_var___u,:internal_var___p,:t]) ==
  """
  void diffeqf(double* internal_var___du, double* internal_var___u, double* internal_var___p, double t) {
    internal_var___du[0] = internal_var___p[0] * internal_var___u[0] + -1 * internal_var___u[0] * internal_var___u[1];
    internal_var___du[1] = internal_var___u[0] * internal_var___u[1] + -3 * internal_var___u[1];
  }
  """

@test Symbolics.build_function(expr,[x,y],[a],t,target = Symbolics.MATLABTarget()) ==
  """
  diffeqf = @(t,internal_var___u) [
    internal_var___p(1) * internal_var___u(1) + -1 * internal_var___u(1) * internal_var___u(2);
    internal_var___u(1) * internal_var___u(2) + -3 * internal_var___u(2);
  ];
  """

@test Symbolics.build_function(expr,[x,y],[a],t,target = Symbolics.CTarget(),
                                     lhsname=:internal_var___du,
                                     rhsnames=[:internal_var___u,:internal_var___p,:t]) ==
    Symbolics.build_function(expr,[x,y],[a],t,target = Symbolics.CTarget(),
                                   lhsname=:internal_var___du,
                                   rhsnames=[:internal_var___u,:internal_var___p,:t])

@test Symbolics.build_function(expr,[x,y],[a],t,target = Symbolics.StanTarget()) ==
      Symbolics.build_function(expr,[x,y],[a],t,target = Symbolics.StanTarget())

@test Symbolics.build_function(expr,[x,y],[a],t,target = Symbolics.MATLABTarget()) ==
      Symbolics.build_function(expr,[x,y],[a],t,target = Symbolics.MATLABTarget())

@test Symbolics.build_function(expr,[x,y],[a],t,target = Symbolics.CTarget(),
                                     lhsname=:internal_var___du,
                                     rhsnames=[:internal_var___u,:internal_var___p,:t]) ==
      Symbolics.build_function(expr,vcat(x,y),[a],t,target = Symbolics.CTarget(),
                                     lhsname=:internal_var___du,
                                     rhsnames=[:internal_var___u,:internal_var___p,:t])

@test Symbolics.build_function(expr,[x,y],[a],t,target = Symbolics.StanTarget()) ==
      Symbolics.build_function(expr,vcat(x,y),[a],t,target = Symbolics.StanTarget())

@test Symbolics.build_function(expr,[x,y],[a],t,target = Symbolics.MATLABTarget()) ==
      Symbolics.build_function(expr,vcat(x,y),[a],t,target = Symbolics.MATLABTarget())


# Matrix CTarget test
let
    @variables x[1:4] y[1:4] z[1:4]
    expression = hcat(x,y,z) # Results in [x1 y1 z1; x2 y2 z2; ...]
    variables  = vcat(x,y,z) # Results in [x1, x2, x3, x4, y1, ...]
    cfunc      = build_function(expression, variables; target = Symbolics.CTarget(), expression = Val{true})

    # Generated function should be out[0] = in[0], out[1] = in[1], out[2] = in[2], etc.
    @test cfunc == "void diffeqf(double* du, double* RHS1) {\n  du[0] = RHS1[0];\n  du[1] = RHS1[1];\n  du[2] = RHS1[2];\n  du[3] = RHS1[3];\n  du[4] = RHS1[4];\n  du[5] = RHS1[5];\n  du[6] = RHS1[6];\n  du[7] = RHS1[7];\n  du[8] = RHS1[8];\n  du[9] = RHS1[9];\n  du[10] = RHS1[10];\n  du[11] = RHS1[11];\n}\n"
end

# Scalar CTarget test
let
    @variables x y z
    expression = x + y + z
    cfunc = build_function(expression, [x], [y], [z]; target = Symbolics.CTarget(), expression = Val{true})

    # Generated function should be out[0] = in1[0] + in2[0] + in3[0]
    @test cfunc == "void diffeqf(double* du, double* RHS1, double* RHS2, double* RHS3) {\n  du[0] = RHS1[0] + RHS2[0] + RHS3[0];\n}\n"
end

# Matrix StanTarget test
let
    @variables x[1:4] y[1:4] z[1:4]
    expression = hcat(x,y,z) # Results in [x1 y1 z1; x2 y2 z2; ...]
    variables  = vcat(x,y,z) # Results in [x1, x2, x3, x4, y1, ...]
    sfunc      = build_function(expression, vcat(x,y), z, []; target = Symbolics.StanTarget(), expression = Val{true})

    # Generated function should be out[1] = v[1], out[2] = v[2], ..., out[9] = p[9], etc.
    @test sfunc == """
    real[] diffeqf(real Any[],real[] internal_var___u,real[] internal_var___p,real[] x_r,int[] x_i) {
      real internal_var___du[12];
      internal_var___du[1] = internal_var___u[1];
      internal_var___du[2] = internal_var___u[2];
      internal_var___du[3] = internal_var___u[3];
      internal_var___du[4] = internal_var___u[4];
      internal_var___du[5] = internal_var___u[5];
      internal_var___du[6] = internal_var___u[6];
      internal_var___du[7] = internal_var___u[7];
      internal_var___du[8] = internal_var___u[8];
      internal_var___du[9] = internal_var___p[1];
      internal_var___du[10] = internal_var___p[2];
      internal_var___du[11] = internal_var___p[3];
      internal_var___du[12] = internal_var___p[4];
      return internal_var___du;
    }
    """
end

# Scalar StanTarget test
let
    @variables t x(t) y(t) z(t)
    expression = x + y + z
    sfunc = build_function(expression, vcat(x,y), [z], t; target = Symbolics.StanTarget(), expression = Val{true})

    # Generated function should be out[0] = v[0] + p[0] + t[0]
    @test sfunc == """
    real[] diffeqf(real t,real[] internal_var___u,real[] internal_var___p,real[] x_r,int[] x_i) {
      real internal_var___du[1];
      internal_var___du[1] = internal_var___u[1] + internal_var___u[2] + internal_var___p[1];
      return internal_var___du;
    }
    """
end

# Matrix MATLABTarget test
let
    @variables x[1:4] y[1:4] z[1:4]
    expression = hcat(x,y,z) # Results in [x1 y1 z1; x2 y2 z2; ...]
    variables  = vcat(x,y,z) # Results in [x1, x2, x3, x4, y1, ...]
    mfunc      = build_function(expression, vcat(x,y,z); target = Symbolics.MATLABTarget(), expression = Val{true})

    # Generated function should be of the same form as expression
    @test mfunc == "diffeqf = @(t,internal_var___u) [\n  internal_var___u(1), internal_var___u(5), internal_var___u(9);\n  internal_var___u(2), internal_var___u(6), internal_var___u(10);\n  internal_var___u(3), internal_var___u(7), internal_var___u(11);\n  internal_var___u(4), internal_var___u(8), internal_var___u(12);\n];\n"
end

# Scalar MATLABTarget test
let
    @variables x y z
    expression = x + y + z
    mfunc = build_function(expression, vcat(x,y,z); target = Symbolics.MATLABTarget(), expression = Val{true})

    # Generated function should be out[0] = v[0] + p[0] + t[0]
    @test mfunc == "diffeqf = @(t,internal_var___u) [\n  internal_var___u(1) + internal_var___u(2) + internal_var___u(3);\n];\n"
end
