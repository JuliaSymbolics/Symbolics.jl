using ModelingToolkit, ParameterizedFunctions

# Build the expected MATLAB string from the actual orderings reported by the
# `System` (parameter index assignment is determined by MTK and changes from
# time to time — e.g. PR #1863's `no-sort-addmul` flipped `[d, c]` to
# `[c, d]`). We look up the index of each named parameter so the test stays
# stable across reorderings.
function _matlab_lotka_volterra_expected(sys, a, b, c, d)
    ps = parameters(sys)
    idx(p) = findfirst(q -> isequal(q, Symbolics.unwrap(p)), Symbolics.unwrap.(ps))
    ia, ib, ic, id = idx(a), idx(b), idx(c), idx(d)
    return "diffeqf = @(internal_var___t,internal_var___u) [\n" *
        "  internal_var___p($ia) * internal_var___u(1) + -1 * internal_var___p($ib) * internal_var___u(1) * internal_var___u(2);\n" *
        "  -1 * internal_var___p($ic) * internal_var___u(2) + internal_var___p($id) * internal_var___u(1) * internal_var___u(2);\n" *
        "];\n"
end

# The order of operands within commutative `*` (and `+`) follows SymbolicUtils'
# argument sort, which is benign but version-dependent (e.g. `u(1) * u(2)` vs
# `u(2) * u(1)`). Compare the two RHS rows modulo that ordering so the test is
# stable across the CI Julia matrix, mirroring the build_targets tests.
function _normalize_matlab(matstr)
    rows = [m[1] for m in eachmatch(r"([^\[\];\n]+);", matstr)]
    return map(rows) do row
        summands = sort(map(strip, split(row, " + ")))
        join(map(s -> join(sort(strip.(split(s, " * "))), " * "), summands), " + ")
    end
end

### Capture MATLABDiffEq.jl type issues
@independent_variables t
@parameters a b c d
@variables x(t) y(t)
D = Differential(t)
eqs = [D(x) ~ a*x - b*x*y
       D(y) ~ -c*y + d*x*y]
@named sys = System(eqs, t)
equations(sys)
matstr = Symbolics.build_function(map(x->x.rhs,equations(sys)),unknowns(sys),
                                        parameters(sys),ModelingToolkit.get_iv(sys),
                                        target = ModelingToolkit.MATLABTarget())
@test _normalize_matlab(matstr) == _normalize_matlab(_matlab_lotka_volterra_expected(sys, a, b, c, d))

f = @ode_def_bare LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c d
p = [1.5,1,3,1]
tspan = (0.0,10.0)
u0 = [1.0,1.0]
prob = ODEProblem(f,u0,tspan,p)

sys = modelingtoolkitize(prob)
[eq.rhs for eq ∈ equations(sys)]
matstr = Symbolics.build_function(map(x->x.rhs,equations(sys)),unknowns(sys),
                                        parameters(sys),ModelingToolkit.get_iv(sys),
                                        target = ModelingToolkit.MATLABTarget())
# `modelingtoolkitize` names the parameters `α, β, γ, δ` (in declaration
# order), so look them up by name from the resulting system.
let ps = parameters(sys)
    a2, b2, c2, d2 = ps[1], ps[2], ps[3], ps[4]
    @test _normalize_matlab(matstr) == _normalize_matlab(_matlab_lotka_volterra_expected(sys, a2, b2, c2, d2))
end
