using Symbolics, Test
using ReferenceTests

# The ordering of operands within commutative `+` and `*` in the generated code
# follows SymbolicUtils' argument sorting, which is benign but has changed across
# versions (e.g. `a * b` vs `b * a`). To keep these tests stable across the CI
# Julia/SymbolicUtils matrix we compare structure modulo that ordering rather than
# byte-for-byte against a reference file. This mirrors the scalar Stan/MATLAB
# tests further down, which were already converted for the same reason.
function normalize_terms(rhs)
    summands = map(strip, split(rhs, " + "))
    canonical = map(summands) do summand
        join(sort(map(strip, split(summand, " * "))), " * ")
    end
    return sort(canonical)
end

@variables t a x(t) y(t)
expr = [a*x - x*y,-3y + x*y]

# StanTarget structural test (was target_functions/1.stan)
let
    sfunc = Symbolics.build_function(expr,[x,y],[a],t,target = Symbolics.StanTarget())
    @test occursin("vector diffeqf(real t,vector internal_var___u,vector internal_var___p)", sfunc)
    @test occursin("vector[2] internal_var___du;", sfunc)
    @test occursin("return internal_var___du;", sfunc)
    m1 = match(r"internal_var___du\[1\] = (.+);", sfunc)
    m2 = match(r"internal_var___du\[2\] = (.+);", sfunc)
    @test m1 !== nothing && m2 !== nothing
    @test normalize_terms(m1[1]) ==
        normalize_terms("internal_var___p[1] * internal_var___u[1] + -1 * internal_var___u[1] * internal_var___u[2]")
    @test normalize_terms(m2[1]) ==
        normalize_terms("-3 * internal_var___u[2] + internal_var___u[1] * internal_var___u[2]")
end

# CTarget structural test (was target_functions/1.c)
let
    cfunc = Symbolics.build_function(expr,[x,y],[a],t,target = Symbolics.CTarget(),
                                     lhsname=:internal_var___du,
                                     rhsnames=[:internal_var___u,:internal_var___p,:t])
    @test occursin("#include <math.h>", cfunc)
    @test occursin("void diffeqf(double* internal_var___du, const double* internal_var___u, const double* internal_var___p, const double t)", cfunc)
    m1 = match(r"internal_var___du\[0\] = (.+);", cfunc)
    m2 = match(r"internal_var___du\[1\] = (.+);", cfunc)
    @test m1 !== nothing && m2 !== nothing
    @test normalize_terms(m1[1]) ==
        normalize_terms("internal_var___p[0] * internal_var___u[0] + -1 * internal_var___u[0] * internal_var___u[1]")
    @test normalize_terms(m2[1]) ==
        normalize_terms("-3 * internal_var___u[1] + internal_var___u[0] * internal_var___u[1]")
end

# MATLABTarget structural test (was target_functions/1.m)
let
    mfunc = Symbolics.build_function(expr,[x,y],[a],t,target = Symbolics.MATLABTarget())
    @test occursin("diffeqf = @(internal_var___t,internal_var___u)", mfunc)
    rows = [m[1] for m in eachmatch(r"([^\[\];\n]+);", mfunc)]
    @test length(rows) == 2
    @test normalize_terms(rows[1]) ==
        normalize_terms("internal_var___p(1) * internal_var___u(1) + -1 * internal_var___u(1) * internal_var___u(2)")
    @test normalize_terms(rows[2]) ==
        normalize_terms("-3 * internal_var___u(2) + internal_var___u(1) * internal_var___u(2)")
end

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
    @test_reference "target_functions/matrix.c" cfunc
end

# Scalar CTarget test
let
    @variables x y t
    expression = x + y + t
    cfunc = build_function(expression, [x], [y], t; target = Symbolics.CTarget(), expression = Val{true})

    @test occursin("void diffeqf(double* du, const double* RHS1, const double* RHS2, const double RHS3)", cfunc)
    m = match(r"du\[0\] = (.+);", cfunc)
    @test m !== nothing
    terms = sort(strip.(split(m[1], "+")))
    @test terms == ["RHS1[0]", "RHS2[0]", "RHS3"]
end

# Scalar CTarget test with scalar multiplication and powers
let
    @variables x y a t
    expression = x^2 + y^-1 + sin(a)^3.5 + 2t + 1//1
    cfunc = build_function(expression, [x, y], [a], t; target = Symbolics.CTarget(), expression = Val{true})

    @test_reference "target_functions/scalar2.c" cfunc
end


# Safe solver functions (ssqrt/scbrt/slog) must map to their math.h names in C
# (https://github.com/JuliaSymbolics/Symbolics.jl/issues/1873)
let
    @variables y
    expression = (1//2) * Symbolics.term(Symbolics.ssqrt, 4y) +
                 Symbolics.term(Symbolics.scbrt, y) + Symbolics.term(Symbolics.slog, y)
    cfunc = build_function([expression], y; target = Symbolics.CTarget(), expression = Val{true})

    @test occursin("sqrt(", cfunc)
    @test occursin("cbrt(", cfunc)
    @test occursin("log(", cfunc)
    @test !occursin("ssqrt", cfunc)
    @test !occursin("scbrt", cfunc)
    @test !occursin("slog", cfunc)
end

# Matrix StanTarget test
let
    @variables x[1:4] y[1:4] z[1:4]
    expression = hcat(x,y,z) # Results in [x1 y1 z1; x2 y2 z2; ...]
    variables  = vcat(x,y,z) # Results in [x1, x2, x3, x4, y1, ...]
    sfunc      = build_function(expression, vcat(x,y), z, []; target = Symbolics.StanTarget(), expression = Val{true})

    @test_reference "target_functions/matrix.stan" sfunc
end

# Scalar StanTarget test
let
    @variables t x(t) y(t) z(t)
    expression = x + y + z
    sfunc = build_function(expression, vcat(x,y), [z], t; target = Symbolics.StanTarget(), expression = Val{true})

    # Term ordering in additions depends on hash-based sorting which is
    # non-deterministic across Julia processes, so we check structurally
    # rather than using an exact reference file.
    @test occursin("vector diffeqf(real t,vector internal_var___u,vector internal_var___p)", sfunc)
    @test occursin("vector[1] internal_var___du;", sfunc)
    @test occursin("return internal_var___du;", sfunc)
    # Check that the RHS contains exactly the expected terms
    m = match(r"internal_var___du\[1\] = (.+);", sfunc)
    @test m !== nothing
    terms = sort(strip.(split(m[1], "+")))
    @test terms == ["internal_var___p[1]", "internal_var___u[1]", "internal_var___u[2]"]
end

# Matrix MATLABTarget test
let
    @variables x[1:4] y[1:4] z[1:4]
    expression = hcat(x,y,z) # Results in [x1 y1 z1; x2 y2 z2; ...]
    variables  = vcat(x,y,z) # Results in [x1, x2, x3, x4, y1, ...]
    mfunc      = build_function(expression, vcat(x,y,z); target = Symbolics.MATLABTarget(), expression = Val{true})

    @test_reference "target_functions/matrix.m" mfunc
end

# Scalar MATLABTarget test
let
    @variables x y z
    expression = x + y + z
    mfunc = build_function(expression, vcat(x,y,z); target = Symbolics.MATLABTarget(), expression = Val{true})

    @test occursin("diffeqf = @(internal_var___t,internal_var___u)", mfunc)
    m = match(r"\[\s*(.+);", mfunc)
    @test m !== nothing
    terms = sort(strip.(split(m[1], "+")))
    @test terms == ["internal_var___u(1)", "internal_var___u(2)", "internal_var___u(3)"]
end
