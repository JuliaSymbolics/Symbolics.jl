"""
    ifelse_eager(cond, x, y)

Symbolic conditional that **always lowers to an eager `ifelse` call** during code
generation. Both branches `x` and `y` are evaluated unconditionally and the result of `cond`
selects between the two computed values, exactly like `Base.ifelse`.

This is the right choice when both branches are always valid to evaluate and you want
branch-free code (e.g. to keep generated code vectorizable/SIMD-friendly, or to avoid
divergence on a GPU). It is contrasted with [`ifelse_branching`](@ref), which emits an
`if`/`else` so that the untaken branch is never evaluated.

`ifelse`, `ifelse_eager` and `ifelse_branching` build the same kind of conditional symbolic
expression and only differ in how [`build_function`](@ref) lowers them. `ifelse` is the
default and is intended to eventually use a cost heuristic to choose a lowering strategy;
the other two pin a specific strategy.

# Examples
```julia
@variables x
f = ifelse_eager(x > 0, x^2, 1 / x)   # both `x^2` and `1/x` are computed, then selected
```

See also [`ifelse_branching`](@ref) and `ifelse`.
"""
function ifelse_eager end

"""
    ifelse_branching(cond, x, y)

Symbolic conditional that **always lowers to an `if`/`else` branch** during code generation,
so that only the taken branch is evaluated. The untaken branch's code is emitted inside the
branch and is therefore never executed when the condition selects the other branch.

This is the right choice when one branch is only valid to evaluate when the condition holds
— for example when a branch divides by a quantity that is zero in the other case, indexes
into something that does not exist, or otherwise errors/produces `NaN`/`Inf` unless its
guard is true. Common subexpression elimination is prevented from hoisting the branches out
of the conditional (see [`build_function`](@ref)'s `cse` option), so the laziness holds even
under CSE. It is contrasted with [`ifelse_eager`](@ref), which evaluates both branches.

`ifelse`, `ifelse_eager` and `ifelse_branching` build the same kind of conditional symbolic
expression and only differ in how [`build_function`](@ref) lowers them. `ifelse` is the
default and is intended to eventually use a cost heuristic to choose a lowering strategy;
the other two pin a specific strategy.

!!! note
    To keep the untaken branch from being evaluated, `ifelse_branching` opts out of common
    subexpression elimination for the whole conditional. As a result, when an
    `ifelse_branching` expression is itself referenced multiple times the generated `if`/`else`
    is duplicated at each use site rather than computed once into a shared temporary (the
    branches remain lazy and the computed values are unchanged). Sharing the conditional
    across use sites while keeping the branches un-hoisted requires a code-generation
    enhancement in SymbolicUtils.

# Examples
```julia
@variables x
# `log(x)` is only evaluated when `x > 0`; the `1/x` branch is only evaluated otherwise.
f = ifelse_branching(x > 0, log(x), 1 / x)
g = build_function(f, x; expression = Val{false})
g(2.0)   # evaluates only `log(2.0)`
```

See also [`ifelse_eager`](@ref) and `ifelse`.
"""
function ifelse_branching end

# Both variants reduce to a plain conditional on concrete (non-symbolic) arguments. By the
# time these methods run the arguments are already evaluated, so the runtime behaviour is
# identical to `ifelse`; the eager/branching distinction only affects generated code.
ifelse_eager(cond, x, y) = Base.ifelse(cond, x, y)
ifelse_branching(cond, x, y) = cond ? x : y

# Type/shape promotion mirrors `ifelse`: the condition must be a `Bool`, both branches must
# share a shape, and the result type is the promotion of the branch types.
function SU.promote_symtype(::typeof(ifelse_eager), C::SU.TypeT, T::SU.TypeT, S::SU.TypeT)
    SU.promote_symtype(ifelse, C, T, S)
end
function SU.promote_symtype(
        ::typeof(ifelse_branching), C::SU.TypeT, T::SU.TypeT, S::SU.TypeT)
    SU.promote_symtype(ifelse, C, T, S)
end
function SU.promote_shape(
        ::typeof(ifelse_eager), shc::SU.ShapeT, sht::SU.ShapeT, shf::SU.ShapeT)
    @nospecialize shc sht shf
    return SU.promote_shape(ifelse, shc, sht, shf)
end
function SU.promote_shape(
        ::typeof(ifelse_branching), shc::SU.ShapeT, sht::SU.ShapeT, shf::SU.ShapeT)
    @nospecialize shc sht shf
    return SU.promote_shape(ifelse, shc, sht, shf)
end

function _conditional_term(op, _if::BasicSymbolic{T}, _then, _else) where {T}
    type = promote_symtype(op, symtype(_if), symtype(_then), symtype(_else))
    sh = SU.promote_shape(op, shape(_if), shape(_then), shape(_else))
    return Term{T}(op, ArgsT{T}((_if, _then, _else)); type, shape = sh)
end

for op in (:ifelse_eager, :ifelse_branching)
    @eval begin
        $op(c::BasicSymbolic, x, y) = _conditional_term($op, c, x, y)

        $op(c::Num, x, y) = $op(unwrap(c), unwrap(x), unwrap(y))
        $op(c::Num, x::Num, y) = Num($op(unwrap(c), unwrap(x), unwrap(y)))
        $op(c::Num, x, y::Num) = Num($op(unwrap(c), unwrap(x), unwrap(y)))
        $op(c::Num, x::Num, y::Num) = Num($op(unwrap(c), unwrap(x), unwrap(y)))

        $op(c::Num, x::Arr{T, N}, y) where {T, N} = Arr{T, N}($op(
            unwrap(c), unwrap(x), unwrap(y)))
        $op(c::Num, x, y::Arr{T, N}) where {T, N} = Arr{T, N}($op(
            unwrap(c), unwrap(x), unwrap(y)))
        $op(c::Num, x::Arr{T, N}, y::Arr{T, N}) where {T, N} = Arr{T, N}($op(
            unwrap(c), unwrap(x), unwrap(y)))
        # Disambiguate mixed scalar/array branches (intersections of the methods above).
        $op(c::Num, x::Num, y::Arr{T, N}) where {T, N} = Arr{T, N}($op(
            unwrap(c), unwrap(x), unwrap(y)))
        $op(c::Num, x::Arr{T, N}, y::Num) where {T, N} = Arr{T, N}($op(
            unwrap(c), unwrap(x), unwrap(y)))
    end
end

# ---- Code generation ----
#
# `ifelse_eager` needs no custom lowering: it is a registered function whose definition calls
# `ifelse`, so the generic fallback lowers it to a plain call that evaluates both branches.
#
# `ifelse_branching` lowers to an `if`/`else` whose branch bodies are generated inside their
# respective arms (via a fresh codegen scope) rather than hoisted before the conditional.
# Combined with the `cse_inside_expr` definition below, this guarantees the untaken branch is
# never evaluated even when CSE is enabled.
function SU.Code.codegen_function!(
        ::typeof(ifelse_branching), cs::SU.Code.CodegenState{T},
        expr::BasicSymbolic{T}, expr_idx::Integer
) where {T}
    cond_idx, true_idx, false_idx = SU.Graphs.outneighbors(cs.ir.dependency_graph, expr_idx)
    cond = cs(cs.ir[cond_idx])

    true_cs, true_bm = SU.Code.enter_scope(cs)
    true_val = true_cs(cs.ir[true_idx])
    true_block = SU.Code.exit_scope!(true_cs, true_bm)
    push!(true_block.args, true_val)

    false_cs, false_bm = SU.Code.enter_scope(cs)
    false_val = false_cs(cs.ir[false_idx])
    false_block = SU.Code.exit_scope!(false_cs, false_bm)
    push!(false_block.args, false_val)

    return SU.Code.codegen!(cs, expr_idx, Expr(:if, cond, true_block, false_block))
end

# Slow `toexpr` path (used by non-`fast_toexpr` build targets). `toexpr` inlines recursively,
# so the branching form is naturally lazy here.
function SU.Code.function_to_expr(::typeof(ifelse_branching), O, st)
    args = arguments(O)
    return Expr(:if, toexpr(args[1], st), toexpr(args[2], st), toexpr(args[3], st))
end

# Prevent CSE from hoisting the branches of `ifelse_branching` out of the conditional, which
# would defeat the laziness. The condition and branches are emitted by the codegen method
# above instead.
SU.Code.cse_inside_expr(sym, ::typeof(ifelse_branching)) = false
