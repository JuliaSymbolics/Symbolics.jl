using Symbolics
using SymbolicUtils
using Test
using Symbolics: value, unwrap, dstar_derivative, dstar_jacobian, jacobian

# Diagnostics: dstar emits factored path sums whose canonical form can differ
# across environments while remaining mathematically equal (e.g. a
# Mul(-1, Add(-x,-y)) vs a flattened Add(x+y)). The assertions below keep the
# original isequal(expand.(...)) semantics, but on failure dump structural
# detail (storage variant, coeff, dict iteration order, cached id/hash fields)
# so CI logs show where the representations diverge.
function _dstar_expr_tree(io, ex, indent)
    ex = Symbolics.unwrap(ex)
    pad = "  "^indent
    if !(ex isa SymbolicUtils.BasicSymbolic)
        println(io, pad, "plain ", typeof(ex), ": ", repr(ex))
        return
    end
    d = getfield(ex, :data)
    variant = split(string(typeof(d)), ".")[end]
    if !SymbolicUtils.iscall(ex)
        println(io, pad, variant, " ", repr(ex))
        return
    end
    print(io, pad, variant, " op=", SymbolicUtils.operation(ex))
    hasproperty(d, :variant) && print(io, " variant=", d.variant)
    hasproperty(d, :coeff) && print(io, " coeff=", repr(d.coeff))
    hasproperty(d, :id) && print(io, " id=", d.id)
    println(io)
    if hasproperty(d, :num)  # Div
        println(io, pad, "  num:")
        _dstar_expr_tree(io, d.num, indent + 2)
        println(io, pad, "  den:")
        _dstar_expr_tree(io, d.den, indent + 2)
    else
        if hasproperty(d, :dict)  # AddMul: expose iteration order + coeffs
            println(io, pad, "  dict order: ",
                join((repr(k) * " => " * repr(v) for (k, v) in d.dict), ", "))
        end
        for c in SymbolicUtils.arguments(ex)
            _dstar_expr_tree(io, c, indent + 1)
        end
    end
end

# Descend two expression trees pairwise and return the first node pair that is
# not isequal while all of its children are — i.e. the exact node whose
# difference makes the whole comparison fail. For AddMul nodes children are
# matched by dict key (order-insensitive) rather than position.
function _first_unequal(a, b)
    a = Symbolics.unwrap(a)
    b = Symbolics.unwrap(b)
    isequal(a, b) && return nothing
    if a isa SymbolicUtils.BasicSymbolic && b isa SymbolicUtils.BasicSymbolic
        da = getfield(a, :data)
        db = getfield(b, :data)
        if typeof(da) === typeof(db)
            if hasproperty(da, :num)  # Div
                r = _first_unequal(da.num, db.num)
                r === nothing && (r = _first_unequal(da.den, db.den))
                r !== nothing && return r
            elseif hasproperty(da, :dict)  # AddMul
                kb = collect(keys(db.dict))
                vb = collect(values(db.dict))
                used = falses(length(kb))
                for (ka, ca) in da.dict
                    idx = findfirst(i -> !used[i] && isequal(kb[i], ka), eachindex(kb))
                    if idx === nothing
                        # no isequal key — descend into the best-matching jac key
                        # (same storage variant; for Divs prefer same denominator)
                        # to expose the inner divergence
                        ua = Symbolics.unwrap(ka)
                        uad = ua isa SymbolicUtils.BasicSymbolic ? getfield(ua, :data) : nothing
                        best = nothing
                        for i in eachindex(kb)
                            used[i] && continue
                            ub = Symbolics.unwrap(kb[i])
                            ub isa SymbolicUtils.BasicSymbolic || continue
                            dbk = getfield(ub, :data)
                            uad !== nothing && typeof(dbk) === typeof(uad) || continue
                            best === nothing && (best = kb[i])
                            # prefer a Div whose denominator matches
                            if hasproperty(dbk, :num) && isequal(dbk.den, uad.den)
                                best = kb[i]
                                break
                            end
                        end
                        best === nothing && return (ka, "no matching key in jac")
                        r = _first_unequal(ka, best)
                        return r === nothing ? (ka, best) : r
                    end
                    used[idx] = true
                    ca == vb[idx] || return (ka, "coeff $ca vs $(vb[idx])")
                    r = _first_unequal(ka, kb[idx])
                    r !== nothing && return r
                end
                da.coeff == db.coeff || return (a, b)
            elseif iscall(a) && iscall(b)  # Term and other arg-based variants
                aa = SymbolicUtils.arguments(a)
                bb = SymbolicUtils.arguments(b)
                if length(aa) == length(bb)
                    for (x, y) in zip(aa, bb)
                        r = _first_unequal(x, y)
                        r !== nothing && return r
                    end
                end
            end
        end
    end
    return (a, b)
end

function _dump_node_fields(io, label, n)
    n = Symbolics.unwrap(n)
    if n isa SymbolicUtils.BasicSymbolic
        d = getfield(n, :data)
        println(io, "  ", label, " variant: ", typeof(d))
        for f in fieldnames(typeof(d))
            f in (:dict, :args, :num, :den) && continue
            v = getfield(d, f)
            f === :id && v === nothing && continue
            println(io, "    ", label, ".", f, " = ", repr(v))
        end
    else
        println(io, "  ", label, ": plain ", typeof(n), " ", repr(n))
    end
end

# dump the post-factoring derivative graph: every node's outgoing edges with
# their values and reachability masks — shows which subgraphs were factored
function _dump_dg(io, roots, vars)
    dg = Symbolics.DerivativeGraph(Symbolics.unwrap.(roots), Symbolics.unwrap.(vars))
    Symbolics.factor_subgraphs!(dg)
    for i in eachindex(dg)
        println(io, "  node ", i, ": ", repr(dg.symbols[i]))
        for e in dg.child_edges[i]
            println(io, "    -> ", e.bott_vertex, "  val=", repr(e.edge_value),
                "  rv=", findall(e.reachable_vars), " rr=", findall(e.reachable_roots))
        end
    end
end

function dstar_eq_jac(dstar_out, jac_out, roots = nothing, vars = nothing)
    d = dstar_out isa AbstractArray ? vec(dstar_out) : [dstar_out]
    j = jac_out isa AbstractArray ? vec(jac_out) : [jac_out]
    ok = true
    for i in eachindex(d, j)
        ed, ej = try
            expand(d[i]), expand(j[i])
        catch err
            println("\n=== expand threw at index ", i, ": ", sprint(showerror, err))
            d[i], j[i]
        end
        isequal(ed, ej) && continue
        ok = false
        println("\n========== dstar/jacobian structural mismatch ==========")
        println("VERSION=", VERSION, " ARCH=", Sys.ARCH, " WORD_SIZE=", Sys.WORD_SIZE)
        println("JULIA_HASH_SEED=", get(ENV, "JULIA_HASH_SEED", "<unset>"))
        println("----- index ", i, " -----")
        println("dstar raw:      ", repr(d[i]))
        println("jac   raw:      ", repr(j[i]))
        println("dstar expanded: ", repr(ed))
        println("jac   expanded: ", repr(ej))
        println("isequal raw:    ", isequal(d[i], j[i]))
        if roots !== nothing
            println("factored derivative graph:")
            try
                _dump_dg(stdout, roots, vars)
            catch err
                println("  dg dump threw: ", sprint(showerror, err))
            end
        end
        try
            println("simplify(expand(d-j)): ", repr(simplify(expand(d[i] - j[i]))))
        catch err
            println("simplify(expand(d-j)) threw: ", sprint(showerror, err))
        end
        println("dstar structure:")
        _dstar_expr_tree(stdout, d[i], 1)
        println("jac structure:")
        _dstar_expr_tree(stdout, j[i], 1)
        println("dstar expanded structure:")
        _dstar_expr_tree(stdout, ed, 1)
        println("jac expanded structure:")
        _dstar_expr_tree(stdout, ej, 1)
        r = _first_unequal(ed, ej)
        if r !== nothing
            x2, y2 = r
            println("first unequal node (expanded): ", repr(x2), "  vs  ", repr(y2))
            _dump_node_fields(stdout, "dstar", x2)
            _dump_node_fields(stdout, "jac", y2)
            r_raw = _first_unequal(d[i], j[i])
            if r_raw !== nothing
                x3, y3 = r_raw
                println("first unequal node (raw): ", repr(x3), "  vs  ", repr(y3))
                _dump_node_fields(stdout, "dstar", x3)
                _dump_node_fields(stdout, "jac", y3)
            end
        end
    end
    return ok
end

# Standard derivative
@variables x
D = Differential(x)

@test isequal(dstar_derivative(x, x), expand_derivatives(D(x)))
@test isequal(dstar_derivative(2x, x), expand_derivatives(D(2x)))
@test isequal(dstar_derivative(x^2, x), expand_derivatives(D(x^2)))
@test isequal(expand(dstar_derivative(sin(cos(x))*cos(cos(x)), x)), expand_derivatives(D(sin(cos(x))*cos(cos(x)))))
@test isequal(expand(dstar_derivative(cos(sin(exp(2x)) + cos(exp(2x))), x)), expand(expand_derivatives(D(cos(sin(exp(2x)) + cos(exp(2x)))))))
@test isequal(expand(dstar_derivative(2sin(x^4 - x) + 3cos(2x), x)), expand(expand_derivatives(D(2sin(x^4 - x) + 3cos(2x)))))
@test dstar_eq_jac(dstar_derivative(2x*exp(x) + 2*exp(x), x), expand(expand_derivatives(D(2x*exp(x) + 2*exp(x)))), [2x*exp(x) + 2*exp(x)], [x])

# Standard Jacobian
@variables x y z
Dx = Differential(x)
Dy = Differential(y)
Dz = Differential(z)

@test isequal(dstar_jacobian([x,y], [x,y]), jacobian([x,y], [x,y]))
@test isequal(dstar_jacobian([x,y,z], [x,y,z]), jacobian([x,y,z], [x,y,z]))
@test isequal(dstar_jacobian([x,y,z], [x]), jacobian([x,y,z], [x]))
@test isequal(dstar_jacobian([x], [x,y,z]), jacobian([x], [x,y,z]))
@test isequal(dstar_jacobian([x*y*z], [x,y,z]), jacobian([x*y*z], [x,y,z]))
@test isequal(dstar_jacobian([x*y*z, x+y+z, sqrt(x^2 + y^2 + z^2)], [x,y,z]), jacobian([x*y*z, x+y+z, sqrt(x^2 + y^2 + z^2)], [x,y,z]))
@test isequal(dstar_jacobian([(x^2+y^2)^2, (x^2+y^2)^2 * y], [x,y]), jacobian([(x^2+y^2)^2, (x^2+y^2)^2 * y], [x,y]))
@test dstar_eq_jac(dstar_jacobian([(x^2 + y^2)*y, (x^2+y^2)*x^2 + (x^2+y^2)*y^2], [x,y]), jacobian([(x^2 + y^2)*y, (x^2+y^2)*x^2 + (x^2+y^2)*y^2], [x,y]), [(x^2 + y^2)*y, (x^2+y^2)*x^2 + (x^2+y^2)*y^2], [x,y])

# Regression tests for parallel edges and edge splitting: factoring creates edges
# that may share endpoints with existing edges and whose reachability extends
# outside the factored subgraph; those outside paths must be preserved
p = x + y
q = x - y
r = p * q
s = p / q
t = r + s
u = r * s
@test dstar_eq_jac(dstar_jacobian([t, u, t * u], [x, y]), jacobian([t, u, t * u], [x, y]), [t, u, t * u], [x, y])
@test dstar_eq_jac(dstar_jacobian([t, u, t * u, t * u + u], [x, y]), jacobian([t, u, t * u, t * u + u], [x, y]), [t, u, t * u, t * u + u], [x, y])
u2 = x^2 + x
@test dstar_eq_jac(dstar_jacobian([u2^2, u2 * x^2], [x]), jacobian([u2^2, u2 * x^2], [x]), [u2^2, u2 * x^2], [x])
@test dstar_eq_jac(dstar_jacobian([u2^2, u2 * x^2, u2 * x^2 + u2], [x]), jacobian([u2^2, u2 * x^2, u2 * x^2 + u2], [x]), [u2^2, u2 * x^2, u2 * x^2 + u2], [x])
# Duplicate arguments: partial derivatives for identical argument positions are
# summed into a single edge
@test isequal(dstar_derivative(atan(x, x), x), expand_derivatives(Differential(x)(atan(x, x))))
@test isequal(dstar_derivative(x^x, x), expand_derivatives(Differential(x)(x^x)))
@test isequal(dstar_derivative(x^x * y + atan(x, x), x), expand_derivatives(Differential(x)(x^x * y + atan(x, x))))

# Edge case Jacobian
@test isequal(dstar_jacobian([x], [x,x]), jacobian([x], [x,x]))
@test isequal(dstar_jacobian([x,x], [x]), jacobian([x,x], [x]))
@test isequal(dstar_jacobian([x,x], [x,x]), jacobian([x,x], [x,x]))

# Duplicate roots and vars are deduplicated before graph construction
@test dstar_eq_jac(dstar_jacobian([t, u, t * u, t * u], [x, y]), jacobian([t, u, t * u, t * u], [x, y]), [t, u, t * u, t * u], [x, y])
@test dstar_eq_jac(dstar_jacobian([u2^2, u2 * x^2, u2^2], [x]), jacobian([u2^2, u2 * x^2, u2^2], [x]), [u2^2, u2 * x^2, u2^2], [x])
@test dstar_eq_jac(dstar_jacobian([t, u], [x, y, x]), jacobian([t, u], [x, y, x]), [t, u], [x, y, x])
@test dstar_eq_jac(dstar_jacobian([t, u, t * u], [x, y, y, x]), jacobian([t, u, t * u], [x, y, y, x]), [t, u, t * u], [x, y, y, x])

# Unregistered functions throw `DerivativeNotDefinedError` instead of asserting
unregistered_fn(a) = a
@register_symbolic unregistered_fn(a)
@test_throws Symbolics.DerivativeNotDefinedError dstar_derivative(unregistered_fn(x), x)

# Array symbolics
@variables z[1:3]
@test isequal(dstar_jacobian(z,z), jacobian(z,z))
@test isequal(dstar_jacobian(2z, z), jacobian(2z, z))

# based on use in FastDifferentiation.jl and the D* paper for testing
function spherical_harmonics(max_l::Integer, x, y, z)
    Pc = Dict{Tuple{Int,Int}, Any}()
    Cc = Dict{Int, Any}()
    Sc = Dict{Int, Any}()

    function P(l, m)
        get!(Pc, (l, m)) do
            if l == 0 && m == 0
                1.0
            elseif l == m
                (1 - 2m) * P(m - 1, m - 1)
            elseif l == m + 1
                (2m + 1) * z * P(m, m)
            else
                ((2l - 1) / (l - m)) * z * P(l - 1, m) - ((l + m - 1) / (l - m)) * P(l - 2, m)
            end
        end
    end

    function C(m)
        get!(Cc, m) do
            m == 0 ? 1 : x * S(m - 1) + y * C(m - 1)
        end
    end

    function S(m)
        get!(Sc, m) do
            m == 0 ? 0 : x * C(m - 1) - y * S(m - 1)
        end
    end

    factorial_approx(n) = sqrt(2π * n) * (n / ℯ * sqrt(n * sinh(1 / n) + 1 / (810 * n^6)))^n
    N(l, m) = m == 0 ? sqrt(2l + 1 / (4π)) : sqrt((2l + 1) / 2π * factorial_approx(l - m) / factorial_approx(l + m))
    Y(l, m) = m < 0 ? N(l, -m) * P(l, -m) * S(-m) : N(l, m) * P(l, m) * C(m)

    return Num[Num(Y(l, m)) for l in 0:max_l-1 for m in -l:l]
end

@variables sx sy sz
sh_vars = [sx, sy, sz]

# max_l=4 -> 16 expressions
sh4 = spherical_harmonics(4, sx, sy, sz)
@test isequal(dstar_jacobian(sh4, sh_vars), jacobian(sh4, sh_vars))

# max_l=5 -> 25 expressions
# (more dominator/postdominator sharing than max_l=4 exercises).
sh5 = spherical_harmonics(5, sx, sy, sz)
@test isequal(dstar_jacobian(sh5, sh_vars), jacobian(sh5, sh_vars))

sh13 = spherical_harmonics(13, sx, sy, sz)
# large enough that dstar_jacobian and jacobian accumulate floating point errors, so sub in vars and use isapprox
let
    dj = dstar_jacobian(sh13, sh_vars)
    j = jacobian(sh13, sh_vars)
    subs = Dict(sx => 0.3, sy => 0.5, sz => 0.7)
    dvals = Symbolics.value.(substitute.(dj, (subs,)))
    jvals = Symbolics.value.(substitute.(j, (subs,)))
    @test isapprox(Float64.(dvals), Float64.(jvals); rtol=1e-8)
end

# Random-DAG fuzzing: random expression trees with shared subexpressions and
# occasional duplicate arguments/roots/vars, compared against jacobian numerically
@testset "fuzz" begin
    using Random
    Random.seed!(7)
    @variables fx fy fz
    fsubs = Dict(fx => 2.3, fy => 0.7, fz => 1.1)

    function rand_expr(depth, pool)
        depth <= 0 && return rand(vcat([fx, fy, fz], pool))
        a = rand_expr(depth - 1, pool)
        b = rand_expr(depth - 1, pool)
        rand() < 0.3 && (b = a) # duplicate arguments exercise identical-child edges
        rand([+, -, *, /])(a, b)
    end
    # expression construction itself can throw (e.g. Num integer-division edge
    # cases); fall back to a variable
    safe_expr(d, p) = try rand_expr(d, p) catch; fx end

    evaluated = 0
    for _ in 1:50
        pool = [safe_expr(rand(1:2), Num[]) for _ in 1:rand(1:3)]
        roots = [safe_expr(rand(2:4), pool) for _ in 1:rand(1:3)]
        rand() < 0.4 && length(roots) > 1 && push!(roots, rand(roots))
        vars = [fx, fy, fz]
        rand() < 0.3 && (vars = [vars; rand(vars)])

        j = try
            jacobian(roots, vars)
        catch
            continue
        end
        dj = dstar_jacobian(roots, vars)
        @test size(dj) == size(j)

        jvals = Float64.(Symbolics.value.(substitute.(j, (fsubs,); fold = Val(true))))
        dvals = Float64.(Symbolics.value.(substitute.(dj, (fsubs,); fold = Val(true))))
        all(isfinite, jvals) && all(isfinite, dvals) || continue
        evaluated += 1
        @test isapprox(dvals, jvals; rtol = 1e-6)
    end
    # guard against the generator producing nothing usable
    @test evaluated > 20
end