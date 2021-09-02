using Symbolics, Test

@variables t a x(t) y(t)
expr = [a*x - x*y,-3y + x*y]
@test Symbolics.build_function(expr,[x,y],[a],t,target = Symbolics.StanTarget()) ==
    """
    vector diffeqf(real t,vector internal_var___u,vector internal_var___p) {
      vector[2] internal_var___du;
      internal_var___du[1] = internal_var___p[1] * internal_var___u[1] + -1 * internal_var___u[1] * internal_var___u[2];
      internal_var___du[2] = internal_var___u[1] * internal_var___u[2] + -3 * internal_var___u[2];
      return internal_var___du;
    }
    """

@test Symbolics.build_function(expr,[x,y],[a],t,target = Symbolics.CTarget(),
                                     lhsname=:internal_var___du,
                                     rhsnames=[:internal_var___u,:internal_var___p,:t]) ==
  """
  #include <math.h>
  void diffeqf(double* internal_var___du, const double* internal_var___u, const double* internal_var___p, const double t) {
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
    @test cfunc == "#include <math.h>\nvoid diffeqf(double* du, const double* RHS1) {\n  du[0] = RHS1[0];\n  du[1] = RHS1[1];\n  du[2] = RHS1[2];\n  du[3] = RHS1[3];\n  du[4] = RHS1[4];\n  du[5] = RHS1[5];\n  du[6] = RHS1[6];\n  du[7] = RHS1[7];\n  du[8] = RHS1[8];\n  du[9] = RHS1[9];\n  du[10] = RHS1[10];\n  du[11] = RHS1[11];\n}\n"
end

# Scalar CTarget test
let
    @variables x y t
    expression = x + y + t
    cfunc = build_function(expression, [x], [y], t; target = Symbolics.CTarget(), expression = Val{true})

    # Generated function should be out[0] = in1[0] + in2[0] + in3
    @test cfunc == "#include <math.h>\nvoid diffeqf(double* du, const double* RHS1, const double* RHS2, const double RHS3) {\n  du[0] = RHS3 + RHS1[0] + RHS2[0];\n}\n"
end

# Scalar CTarget test with scalar multiplication and powers
let
  @variables x y a t
  expression = x^2 + y^-1 + sin(a)^3.5 + 2t
  cfunc = build_function(expression, [x, y], [a], t; target = Symbolics.CTarget(), expression = Val{true})

  # Generated function should avoid scalar multiplication of the form "4t" (currently done by adding another "* 1") and other invalid C syntax
  @test cfunc == "#include <math.h>\nvoid diffeqf(double* du, const double* RHS1, const double* RHS2, const double RHS3) {\n  du[0] = 2 * RHS3 * 1 + pow(RHS1[0], 2) + 1 / RHS1[1] + pow(sin(RHS2[0]), 3.5);\n}\n"
end


# Matrix StanTarget test
let
    @variables x[1:4] y[1:4] z[1:4]
    expression = hcat(x,y,z) # Results in [x1 y1 z1; x2 y2 z2; ...]
    variables  = vcat(x,y,z) # Results in [x1, x2, x3, x4, y1, ...]
    sfunc      = build_function(expression, vcat(x,y), z, []; target = Symbolics.StanTarget(), expression = Val{true})

    # Generated function should be out[1] = v[1], out[2] = v[2], ..., out[9] = p[9], etc.
    @test sfunc == """
    vector diffeqf(real Any[],vector internal_var___u,vector internal_var___p) {
      vector[12] internal_var___du;
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
    vector diffeqf(real t,vector internal_var___u,vector internal_var___p) {
      vector[1] internal_var___du;
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
