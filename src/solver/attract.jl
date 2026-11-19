"""
    detect_addlog(lhs, var)

Traverses expression arg by arg and if two logs of
degree 1 are found it returns true if both of them
have n_func_occ of 1, so they're solvable after 
getting attracted. Returns true if found, 
false if not.

# Arguments
- lhs: A Num expression.
- var: Num and should satisfy is_singleton

# Examples
```jldoctest
julia> using Symbolics

julia> x = (@variables x)[1]
x

julia> Symbolics.detect_addlogs(log(x+1) + log(x-1), x)
true
```
"""
function detect_addlogs(lhs, var)
    u_lhs = unwrap(lhs)
    args = arguments(u_lhs)
    oper = operation(u_lhs)
    !isequal(oper, (+)) && return false

    found = [false, false]
    c = 1
    for arg in args
        !iscall(arg) && continue
        isequal(n_occurrences(arg, var), 0) && continue

        # arg has target var and iscall
        if arg_contains_log(arg, var) && c < 3
            found[c] = true
            c += 1
        end
    end

    return all(found)
end

"""
    detect_exponential(lhs, var)

Traverses an expression and attempts to match the
expression to ``a*b^f(x) + c*d^g(x)``. Hence why it
expects only 2 arguments in the input expression.
Returns true if found, false if not.

# Arguments
- lhs: A Num expression.
- var: Num and should satisfy is_singleton

# Examples
```jldoctest
julia> using Symbolics

julia> x = (@variables x)[1]
x

julia> Symbolics.detect_exponential(2*3^(x+1) + 5^(x^2 + 3), x)
true
```
"""
function detect_exponential(lhs, var)
    lhs = unwrap(lhs)
    !iscall(lhs) && return false
    args = arguments(lhs)
    oper = operation(lhs)
    !isequal(oper, (+)) && return false
    length(args) != 2 && return false

    found = [false, false]
    c = 1
    for arg in args
        n_occurrences(arg, var) == 0 && continue
        !iscall(arg) && continue

        oper = operation(arg)
        if isequal(oper, ^)
            args_arg = arguments(arg)
            n_occurrences(args_arg[1], var) != 0 && continue
            !check_poly_inunivar(args_arg[2], var) && continue

            found[c] = true
            c += 1

        elseif isequal(oper, *) && find_logandexpon(arg, var, ^, 2)
            found[c] = true
            c += 1
        end
    end

    return all(found)
end

"""
    attract_logs(lhs, var)

Rewrites any 2 logs in lhs where the log contains the
input variable var. so ``log(x) + log(y)`` is not attracted.
When attracted, it returns an slog.

# Arguments
- lhs: A Num expression.
- var: Num and should satisfy is_singleton

# Examples
```jldoctest
julia> using Symbolics

julia> x = (@variables x)[1]
x

julia> Symbolics.attract_logs(log(x+1)+log(x-1), x)[1]
slog(-1 + x^2)
```
"""
function attract_logs(lhs, var)
    contains_var(arg) = n_occurrences(arg, var) > 0

    r_addlogs = Vector{Any}()
    r_conditionx = Vector{Any}()
    r_conditiony = Vector{Any}()
    #! format: off
    push!(r_addlogs, @acrule log(~x::(contains_var)) + log(~y::(contains_var)) => slog(~x * ~y))
    push!(r_addlogs, @acrule ~z*log(~x::(contains_var)) + log(~y::(contains_var)) => slog((~x)^(~z) * ~y))
    push!(r_addlogs, @acrule ~z*log(~x::(contains_var)) + ~h*log(~y::(contains_var)) => slog((~x)^(~z) * (~y)^(~h)))

    push!(r_conditionx, @acrule log(~x::(contains_var)) + log(~y::(contains_var)) => ~x)
    push!(r_conditionx, @acrule ~z*log(~x::(contains_var)) + log(~y::(contains_var)) => ~x)
    push!(r_conditionx, @acrule ~z*log(~x::(contains_var)) + ~h*log(~y::(contains_var)) => ~x)
    push!(r_conditiony, @acrule log(~x::(contains_var)) + log(~y::(contains_var)) => ~y)
    push!(r_conditiony, @acrule ~z*log(~x::(contains_var)) + log(~y::(contains_var)) => ~y)
    push!(r_conditiony, @acrule ~z*log(~x::(contains_var)) + ~h*log(~y::(contains_var)) => ~y)
    #! format: on
    condition_x = expand(simplify(lhs, rewriter = SymbolicUtils.Postwalk(SymbolicUtils.Chain(r_conditionx))))
    condition_y = expand(simplify(lhs, rewriter = SymbolicUtils.Postwalk(SymbolicUtils.Chain(r_conditiony))))
    lhs = expand(simplify(lhs, rewriter = SymbolicUtils.Postwalk(SymbolicUtils.Chain(r_addlogs))))

    return lhs, [(condition_x, >), (condition_y, >)]
end

# A key that depends only on the structure of `t`: the arguments of `+` and
# `*` are themselves ordered by their key, so no internal term order leaks in.
# `sorted_arguments` is not enough, because its order ties distinct terms
# such as `2^(x+1)` and `2^(2x-1)` and then keeps the internal order. Leaves
# carry their type and an exact representation, so that values which print
# the same, such as `Float16(0.1)` and `0.1`, still get different keys.
function structural_key(t)
    t = unwrap(t)
    iscall(t) || return leaf_key(SymbolicUtils.unwrap_const(t))
    op = operation(t)
    keys = map(structural_key, arguments(t))
    (op === (+) || op === (*)) && sort!(keys)
    return string(op_key(op), "(", join(keys, ","), ")")
end

# Qualified by module: `Base.log` and `NaNMath.log` both print as `log`. A
# callable object with fields also carries their contents, since two instances
# of one type are different operations.
function op_key(op)
    T = typeof(op)
    op isa Function && Base.issingletontype(T) && return string(parentmodule(op), ".", nameof(op))
    return value_key(op, KEY_DEPTH)
end

# `repr` round-trips the machine floats exactly; a `BigFloat` also needs its
# precision, and a symbol its symbolic type.
leaf_key(v::BasicSymbolic) = string("sym:", v, "::", symtype(v))
leaf_key(v) = value_key(v, KEY_DEPTH)

# How deep `value_key` follows fields, which also stops on cyclic objects.
const KEY_DEPTH = 8

# Key of a non-symbolic value from its contents rather than its display, which
# a custom `show` can make ambiguous: fields are followed recursively, so two
# values that differ in any (nested) field get different keys.
value_key(v::BigFloat, depth) = string("BigFloat@", precision(v), ":", repr(v))
const REPR_EXACT = Union{Base.BitInteger, BigInt, Bool, Base.IEEEFloat, Rational, Complex, AbstractIrrational, AbstractString, Symbol, Char}
value_key(v::REPR_EXACT, depth) = string(typeof(v), ":", repr(v))
function value_key(v::AbstractArray, depth)
    depth == 0 && return string(typeof(v), size(v), "[…]")
    return string(typeof(v), size(v), "[", join((value_key(e, depth - 1) for e in v), ","), "]")
end
function value_key(v, depth)
    T = typeof(v)
    name = string(parentmodule(T), ".", T)
    (nfields(v) == 0 || depth == 0) && return string(name, ":", repr(v))
    fields = (isdefined(v, i) ? value_key(getfield(v, i), depth - 1) : "#undef" for i in 1:nfields(v))
    return string(name, "{", join(fields, ","), "}")
end

# Apply a rule for a two-term sum to the terms of `t` in an order that depends
# only on the terms themselves. `@acrule` tries them in the order of the sum's
# internal term dictionary, which depends on how that sum was first built, so
# the same input could bind `~b` to either base and return a different root.
# `structural_key` gives a readable, stable order, and terms it ranks equally
# are separated by their hash, which is a function of the term's value. Two
# distinct terms that still tie cannot be ordered independently of the
# internal order, so such a sum is left alone rather than given a root picked
# by that order. Since the key follows the contents of every leaf and callable,
# such terms have bases and exponents with equal contents, and for pure
# callables `b^f + b^f = 0` has no solution, so nothing valid is lost. Like
# `@acrule`, the terms the rule does not consume are kept.
function apply_to_sorted_pairs(rule, t)
    (iscall(t) && operation(t) === (+)) || return nothing
    args = collect(arguments(t))
    order_keys = [(structural_key(a), hash(a)) for a in args]
    perm = sortperm(order_keys)
    args, order_keys = args[perm], order_keys[perm]
    n = length(args)
    # The terms of a sum are distinct, so equal keys are an unresolved tie.
    any(i -> order_keys[i] == order_keys[i + 1], 1:(n - 1)) && return nothing
    for i in 1:n, j in 1:n
        i == j && continue
        res = rule(term(+, args[i], args[j]))
        res === nothing && continue
        rest = [args[k] for k in 1:n if k != i && k != j]
        return isempty(rest) ? res : +(res, rest...)
    end
    return nothing
end

"""
    attract_exponential(lhs, var)

Rewrites ``a*b^f(x) + c*d^g(x)`` into 
``f(x) * log(b) - g(x) * log(d) + log(-a/c)``
Which is a solvable polynomial, hence reducing the smart
number of occurrences of x in the expression.

# Arguments
- lhs: A Num expression.
- var: Num and should satisfy is_singleton

# Examples
```jldoctest
julia> using Symbolics

julia> x = (@variables x)[1]
x

julia> Symbolics.attract_exponential(2^(x+1) + 5^(x+3), x)
slog(-1) + slog(2) - 3slog(5) + slog(2)*x - slog(5)*x
```
"""
function attract_exponential(lhs, var)
    lhs = unwrap(lhs)
    contains_var(arg) = n_occurrences(arg, var) > 0

    r_addexpon = Vector{Any}()

    #! format: off
    push!(r_addexpon, @rule (~b)^(~f::(contains_var)) + (~d)^(~g::(contains_var)) => ~f*term(slog, ~b) - ~g*term(slog, ~d) + term(slog, -1))
    push!(r_addexpon, @rule (~a)*(~b)^(~f::(contains_var)) + (~d)^(~g::(contains_var)) => ~f*term(slog, ~b) - ~g*term(slog, ~d) + term(slog, -~a))
    push!(r_addexpon, @rule (~a)*(~b)^(~f::(contains_var)) + (~c)*(~d)^(~g::(contains_var)) => ~f*term(slog, ~b) - ~g*term(slog, ~d) + term(slog, -(~a)//(~c)))
    #! format: on
    r_addexpon = [t -> apply_to_sorted_pairs(r, t) for r in r_addexpon]

    lhs = expand(simplify(
        lhs, rewriter = SymbolicUtils.Postwalk(SymbolicUtils.Chain(r_addexpon))))

    return expand(lhs)
end

"""
    attract_trig(lhs, var)

This function brute forces a bunch of trig identities
on the input lhs as an attempt to reduce the occurrences
of x in the lhs.

# Arguments
- lhs: A Num expression.
- var: Num and should satisfy is_singleton

# Examples
```jldoctest
julia> using Symbolics

julia> x = (@variables x)[1]
x

julia> Symbolics.attract_trig(2sin(x)cos(x), x)
sin(2x)

julia> Symbolics.attract_trig(sin(x)^2 + cos(x)^2, x)
1

julia> Symbolics.attract_trig(cosh(x)^2 + sinh(x)^2, x)
cosh(2x)
"""
function attract_trig(lhs, var)
    lhs = unwrap(lhs)
    contains_var(arg) = n_occurrences(arg, var) > 0

    # r_doubleangle1 = @acrule 2*sin(~x::(contains_var))*cos(~x::(contains_var)) => sin(2*~x)
    r_trig = [@acrule(sin(~x::(contains_var))^2 + cos(~x::(contains_var))^2=>one(~x))
              @acrule(sin(~x::(contains_var))^2 + -1=>-1 * cos(~x)^2)
              @acrule(cos(~x::(contains_var))^2 + -1=>-1 * sin(~x)^2)
              @acrule(cos(~x::(contains_var))^2 + -1 * sin(~x::(contains_var))^2=>cos(2*~x))
              @acrule(sin(~x::(contains_var))^2 + -1 * cos(~x::(contains_var))^2=>-cos(2*~x))
              @acrule(cos(~x::(contains_var)) * sin(~x::(contains_var))=>sin(2 * ~x) / 2)
              @acrule(tan(~x::(contains_var))^2 + -1 * sec(~x::(contains_var))^2=>one(~x))
              @acrule(-1 * tan(~x::(contains_var))^2 + sec(~x::(contains_var))^2=>one(~x))
              @acrule(tan(~x::(contains_var))^2 + 1=>sec(~x)^2)
              @acrule(sec(~x::(contains_var))^2 + -1=>tan(~x)^2)
              @acrule(cot(~x::(contains_var))^2 + -1 * csc(~x)^2=>one(~x))
              @acrule(cot(~x::(contains_var))^2 + 1=>csc(~x)^2)
              @acrule(csc(~x::(contains_var))^2 + -1=>cot(~x)^2)
              @acrule(cosh(~x::(contains_var))^2 + -1 * sinh(~x)^2=>one(~x))
              @acrule(cosh(~x::(contains_var))^2 + -1=>sinh(~x)^2)
              @acrule(sinh(~x::(contains_var))^2 + 1=>cosh(~x)^2)
              @acrule(cosh(~x::(contains_var))^2 + sinh(~x::(contains_var))^2=>cosh(2 * ~x))
              @acrule(cosh(~x::(contains_var)) * sinh(~x::(contains_var))=>sinh(2 * ~x) / 2)]

    lhs = expand(simplify(
        lhs, rewriter = SymbolicUtils.Postwalk(SymbolicUtils.Chain(r_trig))))

    return lhs
end
