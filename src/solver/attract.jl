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
julia> detect_addlogs(log(x+1) + log(x-1), x)
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
julia> detect_exponential(2*3^(x+1) + 5^(x^2 + 3), x)
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
julia> RootFinding.attract_logs(log(x+1)+log(x-1), x)
Main.RootFinding.slog(-1 + x^2)
```
"""
function attract_logs(lhs, var)
    contains_var(arg) = n_occurrences(arg, var) > 0

    r_addlogs = Vector{Any}()
    #! format: off
    push!(r_addlogs, @acrule log(~x::(contains_var)) + log(~y::(contains_var)) => slog(~x * ~y))
    push!(r_addlogs, @acrule ~z*log(~x::(contains_var)) + log(~y::(contains_var)) => slog((~x)^(~z) * ~y))
    push!(r_addlogs, @acrule ~z*log(~x::(contains_var)) + ~h*log(~y::(contains_var)) => slog((~x)^(~z) * (~y)^(~h)))
    #! format: on

    lhs = expand(simplify(
        lhs, rewriter = SymbolicUtils.Postwalk(SymbolicUtils.Chain(r_addlogs))))

    return lhs
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
julia> attract_exponential(2^(x+1) + 5^(x+3), x)
Main.RootFinding.slog(2) + log(complex(-1)) - 3Main.RootFinding.slog(5) + x*Main.RootFinding.slog(2) - x*Main.RootFinding.slog(5)
```
"""
function attract_exponential(lhs, var)
    lhs = unwrap(lhs)
    contains_var(arg) = n_occurrences(arg, var) > 0

    r_addexpon = Vector{Any}()

    #! format: off
    push!(r_addexpon, @acrule (~b)^(~f::(contains_var)) + (~d)^(~g::(contains_var)) => ~f*term(slog, ~b) - ~g*term(slog, ~d) + term(log, term(complex, -1)))
    push!(r_addexpon, @acrule (~a)*(~b)^(~f::(contains_var)) + (~d)^(~g::(contains_var)) => ~f*term(slog, ~b) - ~g*term(slog, ~d) + term(slog, -~a))
    push!(r_addexpon, @acrule (~a)*(~b)^(~f::(contains_var)) + (~c)*(~d)^(~g::(contains_var)) => ~f*term(slog, ~b) - ~g*term(slog, ~d) + term(slog, -(~a)//(~c)))
    #! format: on

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
julia> attract_trig(2sin(x)cos(x), x)
sin(2x)

julia> attract_trig(sin(x)^2 + cos(x)^2, x)
1

julia> RootFinding.attract_trig(cosh(x)^2 + sinh(x)^2, x)
cosh(2x)
"""
function attract_trig(lhs, var)
    lhs = unwrap(lhs)
    contains_var(arg) = n_occurrences(arg, var) > 0

    # r_doubleangle1 = @acrule 2*sin(~x::(contains_var))*cos(~x::(contains_var)) => sin(2*~x)
    r_trig = [@acrule(sin(~x::(contains_var))^2 + cos(~x::(contains_var))^2=>one(~x))
              @acrule(sin(~x::(contains_var))^2 + -1=>-1 * cos(~x)^2)
              @acrule(cos(~x::(contains_var))^2 + -1=>-1 * sin(~x)^2)
              @acrule(cos(~x::(contains_var))^2 + -1 * sin(~x::(contains_var))^2=>cos(2 *
                                                                                      ~x))
              @acrule(sin(~x::(contains_var))^2 + -1 * cos(~x::(contains_var))^2=>-cos(2 *
                                                                                       ~x))
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
