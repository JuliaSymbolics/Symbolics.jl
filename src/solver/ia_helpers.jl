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
    try
        args = arguments(argument)
    catch e
        if isequal(argument, var)
            return 1 
        end
        return 0
    end

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
    counted_ops = [sqrt, cbrt, sin, log, log2, log10, cos, tan, asin, acos, atan, exp, ssqrt, scbrt, slog]
    n = 0


    if cur_oper === (*) || cur_oper === (+)

        outside = false
        for arg in args
            n_occurrences(arg, var) == 0 && continue

            if !iscall(arg) && isequal(var, get_variables(expr)[1]) && !outside
                n += 1
                outside = true
                continue
            end
            !iscall(arg) && continue
            oper = operation(arg)

            args_arg = arguments(arg)
            oper_arg = operation(arg)
            is_var_outside(arg) = check_poly_inunivar(arg, var) && !outside && n_occurrences(arg, var) != 0
            case_1_pow = oper_arg === (^) && n_occurrences(args_arg[2], var) == 0  && n_occurrences(args_arg[1], var) != 0 && is_var_outside(arg)
            case_2_pow = oper_arg === (^) && n_occurrences(args_arg[2], var) != 0  && n_occurrences(args_arg[1], var) == 0  
            case_3_pow = oper_arg === (^) && n_occurrences(args_arg[2], var) == 0  && n_occurrences(args_arg[1], var) != 0 && !check_poly_inunivar(args_arg[1], var)

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
