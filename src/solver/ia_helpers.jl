function n_occurrences(expr, var)
    n = 0
    !iscall(unwrap(expr)) && any(isequal(var, x) for x in get_variables(expr)) && return 1
    !iscall(unwrap(expr)) && return 0

    args = arguments(unwrap(expr))

    for arg in args
        n += traverse(arg, var)
    end

    return n
end

function traverse(argument, var)
    args = []
    !iscall(argument) && isequal(argument, var) && return 1
    !iscall(argument) && !isequal(argument, var) && return 0

    args = arguments(argument)
    n = 0

    for arg in args
        n += traverse(arg, var)
    end
    return n
end

function n_func_occ(expr, var)
    expr = unwrap(expr)
    !iscall(expr) && return n_occurrences(expr, var)
    args, cur_oper = arguments(expr), operation(expr)
    counted_ops = [sqrt, cbrt, sin, log, log2, log10, cos, tan,
        asin, acos, atan, exp, ssqrt, scbrt, slog]
    n = 0

    if cur_oper === (*) || cur_oper === (+)
        outside = false
        for arg in args
            n_occurrences(arg, var) == 0 && continue

            # x
            if !iscall(arg) && isequal(var, get_variables(arg)[1]) && !outside
                n += 1
                outside = true
                continue
            end

            !iscall(arg) && continue
            oper = operation(arg)

            args_arg = arguments(arg)
            oper_arg = operation(arg)
            function is_var_outside(arg)
                check_poly_inunivar(arg, var) && !outside && n_occurrences(arg, var) != 0
            end
            case_1_pow = oper_arg === (^) && n_occurrences(args_arg[2], var) == 0 &&
                         n_occurrences(args_arg[1], var) != 0 &&
                         check_poly_inunivar(args_arg[1], var) &&
                         n_occurrences(arg, var) != 0 && !(args_arg[2] isa Number)
            case_2_pow = oper_arg === (^) && n_occurrences(args_arg[2], var) != 0 &&
                         n_occurrences(args_arg[1], var) == 0
            case_3_pow = oper_arg === (^) && n_occurrences(args_arg[2], var) == 0 &&
                         n_occurrences(args_arg[1], var) != 0 &&
                         !check_poly_inunivar(args_arg[1], var)

            # any transcedental operation and the case:  (weird_transcedental_f(x))^(something)
            if any(isequal(oper, op) for op in counted_ops) || case_3_pow
                n += n_func_occ(args_arg[1], var)

                # the case (some constant)^(f(x))
            elseif case_2_pow
                n += n_func_occ(args_arg[2], var)

                # var is outside 'x'+1
            elseif is_var_outside(arg)
                n += 1
                outside = true

                # case (f(x))^(weird stuff)
            elseif case_1_pow
                n += 1

                # n(2 / x) = 1; n(x/x^2) = 2?
            elseif oper_arg === (/)
                n += n_func_occ(args_arg[1], var)
                n += n_func_occ(args_arg[2], var)

                # multiplication cases
            elseif oper_arg === (*)
                args_arg = arguments(arg)
                for sub_arg in args_arg
                    # x*log(2)
                    if is_var_outside(sub_arg)
                        n += 1
                        outside = true
                        # log(x)*y
                    elseif !check_poly_inunivar(sub_arg, var)
                        n += n_func_occ(sub_arg, var)
                    end
                end
            end

        end

    else
        for arg in args
            n += n_func_occ(arg, var)
        end
    end

    return n
end

function arg_contains_log(arg, var)
    oper = operation(arg)
    isequal(oper, log) && return true

    if oper == (*)
        return find_logandexpon(arg, var, log, 1)
    end

    return false
end

function find_logandexpon(arg, var, oper, poly_index)
    args_arg = arguments(arg)

    oper_term, constant_term = 0, 0

    for a in args_arg
        if n_occurrences(a, var) != 0 && iscall(a) && operation(a) == (oper) &&
           check_poly_inunivar(arguments(a)[poly_index], var)
            oper_term = a
        elseif n_occurrences(a, var) == 0
            constant_term = a
        end
    end

    !isequal(oper_term, 0) && !isequal(constant_term, 0) && return true
    return false
end

"""
    ia_conditions!(f, lhs, rhs::Vector{Any}, conditions::Vector{Tuple})

If `f` is a left-invertible function, `lhs` and `rhs[i]` are univariate functions and
`f(lhs) ~ rhs[i]` for all `i in eachindex(rhss)`, push to `conditions` all the relevant
conditions on `lhs` or `rhs[i]`. Each condition is of the form `(sym, op)` where `sym`
is an expression involving `lhs` and/or `rhs[i]` and `op` is a binary relational operator.
The condition `op(sym, 0)` is then required to be true for the equation `f(lhs) ~ rhs[i]`
to be valid.

For example, if `f = log`, `lhs = x` and `rhss = [y, z]` then the condition `x > 0` must
be true. Thus, `(lhs, >)` is pushed to `conditions`. Similarly, if `f = sqrt`, `rhs[i] >= 0`
must be true for all `i`, and so `(y, >=)` and `(z, >=)` will be appended to `conditions`. 
"""
function ia_conditions!(args...; kwargs...) end

for fn in [log, log2, log10, NaNMath.log, NaNMath.log2, NaNMath.log10, slog]
    @eval function ia_conditions!(::typeof($fn), lhs, rhs, conditions)
        push!(conditions, (lhs, >))
    end
end

for fn in [log1p, NaNMath.log1p]
    @eval function ia_conditions!(::typeof($fn), lhs, rhs, conditions)
        push!(conditions, (lhs - 1, >))
    end
end

for fn in [sqrt, NaNMath.sqrt, ssqrt]
    @eval function ia_conditions!(::typeof($fn), lhs, rhs, conditions)
        for r in rhs
            push!(conditions, (r, >=))
        end
    end
end

"""
    is_periodic(f)

Return `true` if `f` is a single-input single-output periodic function. Return `false` by
default. If `is_periodic(f) == true`, then `fundamental_period(f)` must also be defined.

See also: [`fundamental_period`](@ref)
"""
is_periodic(f) = false

for fn in [
    sin, cos, tan, csc, sec, cot, NaNMath.sin, NaNMath.cos, NaNMath.tan, sind, cosd, tand,
    cscd, secd, cotd, cospi
]
    @eval is_periodic(::typeof($fn)) = true
end

"""
    fundamental_period(f)

Return the fundamental period of periodic function `f`. Must only be called if
`is_periodic(f) == true`.

see also: [`is_periodic`](@ref)
"""
function fundamental_period end

for fn in [sin, cos, csc, sec, NaNMath.sin, NaNMath.cos]
    @eval fundamental_period(::typeof($fn)) = 2pi
end

for fn in [sind, cosd, cscd, secd]
    @eval fundamental_period(::typeof($fn)) = 360.0
end

fundamental_period(::typeof(cospi)) = 2.0

for fn in [tand, cotd]
    @eval fundamental_period(::typeof($fn)) = 180.0
end

for fn in [tan, cot, NaNMath.tan]
    # `1pi isa Float64` whereas `pi isa Irrational{:π}`
    @eval fundamental_period(::typeof($fn)) = 1pi
end
