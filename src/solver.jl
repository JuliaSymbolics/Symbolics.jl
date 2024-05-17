using Symbolics
using LambertW

#  =======  MAIN FUNCTIONS ======

#=
returns solution/s to the equation in terms of the variable

single_solution sets weather it returns only one solution or a set
=#
function solve_single_eq(
    eq::Equation,
    var,
    single_solution = false,
    verify = true,
)
    unchecked_solutions = solve_single_eq_unchecked(eq, var, single_solution)

    unchecked_solutions == nothing && return nothing

    !verify && return unchecked_solutions

    unchecked_solutions =
        !(unchecked_solutions isa Vector) ? [unchecked_solutions] : unchecked_solutions
    try
        float_solutions = convert_solutions_to_floats(unchecked_solutions)

        if (float_solutions != nothing)#check answers to make sure they are valid solutions numerically
            to_be_removed = falses(0)

            for float_solution in float_solutions
                left_side_fval =
                    convert(Float64, substitute(eq.lhs, [var => float_solution]))
                right_side_fval =
                    convert(Float64, substitute(eq.rhs, [var => float_solution]))
                push!(
                    to_be_removed,
                    !(abs(left_side_fval - right_side_fval) <= 1.0 / 2.0^20.0),
                )
            end

            deleteat!(unchecked_solutions, to_be_removed)
        end
    catch
        @warn "unable to verify solutions"
    end

    if length(unchecked_solutions) == 1
        unchecked_solutions = unchecked_solutions[1]
    end

    return unchecked_solutions
end

#returns a dictionary of the solutions to the system of equations
function solve_system_eq(equs::Vector{Equation}, vars)
    removed = Equation[]#keep track of removed equations

    reduced = copy(equs)#the reducing set

    for i = 1:length(vars)#go through each variable and remove
        remove_eq(reduced, vars[i], removed)
    end

    #solve last equation in the reduced set

    solutions = Dict()

    #re substitute the variables back in to find value
    for i = length(removed):-1:1
        current_eq = substitute(removed[i], solutions)
        solutions[current_eq.lhs] = current_eq.rhs
    end

    return solutions
end

#  =======  MAIN FUNCTIONS END ======

function get_parts_list(a, b, a_list = Vector{Any}(), b_list = Vector{Any}())
    if SymbolicUtils.issym(a)
        push!(a_list, a)
        push!(b_list, b)
    elseif iscall(a) && iscall(b) && isequal(operation(a), operation(b))
        a_args = arguments(a)
        b_args = arguments(b)

        length(a_args) != length(b_args) && return Nothing

        for i = 1:length(a_args)
            check = get_parts_list(a_args[i], b_args[i], a_list, b_list)
            check == Nothing && return Nothing
        end
    elseif a isa Equation && b isa Equation
        get_parts_list(a.lhs, b.lhs, a_list, b_list)
        get_parts_list(a.rhs, b.rhs, a_list, b_list)
    else
        return Nothing
    end
    return (a_list, b_list)
end

function find_matches(a)
    vars_seen = UInt[]
    matching = Vector{Int}[]

    for i = 1:length(a)
        match_set = Int[]

        push!(match_set, i)

        current = a[i]
        (hash(current) in vars_seen) && continue

        for j = i+1:length(a)
            other = a[j]
            if isequal(current, other)
                push!(match_set, j)
            end
        end
        push!(vars_seen, hash(current))

        if length(match_set) > 1
            push!(matching, match_set)
        end

    end
    return matching
end

function verify_matches(a_parts, b_parts)
    a_matches = find_matches(a_parts)
    for match in a_matches
        to_compare = b_parts[match[1]]
        for i = 2:length(match)
            !isequal(to_compare, b_parts[match[i]]) && return false
        end
    end
    return true
end

struct Becomes
    before::Any
    after::Any
end

function create_eq_pairs(a, b)
    out = Dict{BasicSymbolic,Any}()

    for i = 1:length(a)
        a_part = a[i]
        a_part in keys(out) && continue

        out[a_part] = b[i]
    end

    return out
end

function replace_term(expr, dic::Dict)
    if SymbolicUtils.issym(expr) && haskey(dic, expr)
        return dic[expr]
    elseif iscall(expr)
        args = Any[]

        for arg in arguments(expr)
            push!(args, replace_term(arg, dic))
        end

        op = operation(expr)

        if SymbolicUtils.isterm(expr)
            return term(op, args...)
        else
            op == (/) && return SymbolicUtils.Div(args[1], args[2])
            return op(args...)
        end
    elseif expr isa Equation
        return replace_term(expr.lhs, dic) ~ replace_term(expr.rhs, dic)
    else
        return expr
    end
end

function apply_ns_rule(becomes::Becomes, expr)
    if expr_similar(becomes.before, expr, false)
        check = get_parts_list(becomes.before, expr)
        check == Nothing && return expr

        becomes_parts, expr_parts = check

        !verify_matches(becomes_parts, expr_parts) && return expr

        eq_list = create_eq_pairs(becomes_parts, expr_parts)

        return replace_term(becomes.after, eq_list)

    end
    return expr
end

function expr_similar(ref_expr, expr, check_matches = true)

    SymbolicUtils.issym(ref_expr) && return true
    SymbolicUtils.issym(expr) && iscall(ref_expr) && return false

    if iscall(ref_expr)
        ref_args = arguments(ref_expr)
        ref_len = length(ref_args)
        ref_op = operation(ref_expr)

        args = arguments(expr)
        len = length(args)
        op = operation(expr)

        (!isequal(ref_op, op) || !isequal(ref_len, len)) && return false

        for i = 1:ref_len
            ref_arg = ref_args[i]
            arg = args[i]

            !expr_similar(ref_arg, arg, false) && return false
        end

        if check_matches
            becomes_parts, expr_parts = get_parts_list(ref_expr, expr)
            !verify_matches(becomes_parts, expr_parts) && return false
        end

        return true
    elseif ref_expr isa Equation && expr isa Equation

        if (check_matches)
            check = get_parts_list(ref_expr, expr)
            check == Nothing && return false
            ref_parts, expr_parts = check

            !verify_matches(ref_parts, expr_parts) && return false
        end

        return expr_similar(ref_expr.lhs, expr.lhs, false) &&
               expr_similar(ref_expr.rhs, expr.rhs, false) &&
               return true
    else
        return isequal(ref_expr, expr)
    end

    return false
end

function get_base(expr)
    (!iscall(expr) || operation(expr) != (^)) && throw("not a power(^) -> $expr")
    return arguments(expr)[1]
end

function get_exp(expr)
    (!iscall(expr) || operation(expr) != (^)) && throw("not a power(^) -> $expr")
    return arguments(expr)[2]
end

function solve_single_eq_unchecked(
    eq::Equation,
    var,
    single_solution = false,
)
    eq = (SymbolicUtils.add_with_div(eq.lhs + -1 * eq.rhs) ~ 0)#move everything to the left side

    #eq = termify(eq)

    while (true)
        oldState = eq

        if (iscall(eq.lhs))

            potential_solution = solve_quadratic(eq, var, single_solution)
            if potential_solution isa Equation
                eq = potential_solution
            else
                return potential_solution
            end


            op = operation(eq.lhs)

            if (op in (+, *))#N argumented types

                eq = move_to_other_side(eq, var)
                eq = special_strategy(eq, var)

                potential_solution = left_prod_right_zero(eq, var, single_solution)
                if potential_solution isa Equation
                    eq = potential_solution
                else
                    return potential_solution
                end

            elseif (op == /)#reverse division
                eq = eq.lhs.num - eq.lhs.den * eq.rhs ~ 0
            elseif (op == ^)#reverse powers

                potential_solution =
                    reverse_powers(eq::Equation, var, single_solution)
                if potential_solution isa Equation
                    eq = potential_solution
                else
                    return potential_solution
                end

            else
                eq = inverse_funcs(eq::Equation, var)
            end

        end
        if (isequal(eq.lhs, var))
            return eq#solved!
        end


        if (isequal(eq, oldState))
            @warn "unable to solve $(eq) in terms of $var"
            return nothing#unsolvable with these methods
        end
    end
end

function left_prod_right_zero(eq::Equation, var, single_solution)
    if SymbolicUtils.ismul(eq.lhs) && isequal(0, eq.rhs)
        if (single_solution)
            eq = arguments(eq.lhs)[1] ~ 0
        else
            solutions = Equation[]
            for arg in arguments(eq.lhs)
                temp = solve_single_eq(arg ~ 0, var)
                temp = temp isa Equation ? [temp] : temp
                push!(solutions, temp...)
            end
            return solutions
        end
    end
    return eq
end

#reduce the system of equations by one equation according to the provided variable
function remove_eq(equs::Vector{Equation}, var, removed::Vector{Equation})
    for i = 1:length(equs)
        solution = solve_single_eq(equs[i], var, true)
        if (solution isa Vector)
            solution = solution[1]
        end
        solution == nothing && continue

        push!(removed, solution)
        deleteat!(equs, i)

        for j = 1:length(equs)
            equs[j] = substitute(equs[j], Dict(solution.lhs => solution.rhs))
        end
        return
    end
end
#=
moves non variable components to the other side of the equation

example move_to_other_side(x+a~z,x) returns x~z-a

=#
function move_to_other_side(eq::Equation, var)

    !iscall(eq.lhs) && return eq#make sure left side is tree form

    op = operation(eq.lhs)

    if op in (+, *)
        elements = arguments(eq.lhs)

        stays = []#has variable
        move = []#does not have variable
        for i = 1:length(elements)
            hasVar = SymbolicUtils._occursin(var, elements[i])
            if (hasVar)
                push!(stays, elements[i])
            else
                push!(move, elements[i])
            end
        end

        if (op == +)#reverse addition
            eq =
                (length(stays) == 0 ? 0 : +(stays...)) ~
                    -(length(move) == 0 ? 0 : +(move...)) + eq.rhs
        elseif (op == *)#reverse multiplication
            eq =
                (length(stays) == 0 ? 1 : *(stays...)) ~
                    SymbolicUtils.Div(eq.rhs, (length(move) == 0 ? 1 : *(move...)))
        end
    end
    return eq
end

#more rare solving strategies
function special_strategy(eq::Equation, var)

    @syms a b c x y z
    rules = [
        Becomes(
            sin(x) + cos(x) ~ y,
            x ~ acos(y / term(sqrt, 2)) + SymbolicUtils.Div(pi, 4),
        ),
        Becomes(
            x * a^x ~ y,
            x ~ SymbolicUtils.Div(term(lambertw, y * term(log, a)), term(log, a)),
        ),
        Becomes(x * log(x) ~ y, x ~ SymbolicUtils.Div(y, term(lambertw, y))),
        Becomes(x * exp(x) ~ y, x ~ term(lambertw, y)),
        Becomes(a + sqrt(b) ~ c, b - a^2 + 2 * a * c - c^2 ~ 0),
    ]


    for rule in rules
        eq = apply_ns_rule(rule, eq)
    end


    !iscall(eq.lhs) && return eq#make sure left side is tree form

    op = operation(eq.lhs)
    elements = arguments(eq.lhs)

    if (op == +) &&
       length(elements) == 2 &&
       sum(iscall.(elements)) == length(elements) &&
       isequal(operation.(elements), [sqrt for el = 1:length(elements)]) #check for sqrt(a)+sqrt(b)=c form , to solve this sqrt(a)+sqrt(b)=c -> 4*a*b-full_expand((c^2-b-a)^2)=0 then solve using quadratics

        #grab values
        a = (elements[1]).arguments[1]
        b = (elements[2]).arguments[1]
        c = eq.rhs


        eq =
            expand(2 * b * a) - expand(a^2) - expand(b^2) - expand(c^4) +
            expand(2 * a * c^2) +
            expand(2 * b * c^2) ~ 0

    elseif (op == +) &&
           isequal(eq.rhs, 0) &&
           length(elements) == 2 &&
           sum(iscall.(elements)) == length(elements) &&
           length(arguments(elements[1])) == 2 &&
           isequal(arguments(elements[1])[1], -1) &&
           iscall(arguments(elements[1])[2]) &&
           operation(elements[2]) == operation(arguments(elements[1])[2])#-f(y)+f(x)=0 -> x-y=0

        x = arguments(elements[2])[1]
        y = arguments(arguments(elements[1])[2])[1]

        eq = x - y ~ 0
    end

    return eq
end

#=
reduces the form of the square root
examples
reduce_root(term(sqrt,64)) = 8
reduce_root(term(sqrt,32)) = 4*sqrt(2)
=#

function reduce_root(a)

    if SymbolicUtils.ispow(a) && a.exp isa Rational
        a = term(^, a.base, a.exp)
    end

    if iscall(a) && (operation(a) == sqrt)
        a = SymbolicUtils.Pow(arguments(a)[1], 1 // 2)
    elseif iscall(a) &&
           (operation(a) == ^) &&
           isequal(arguments(a)[2], 1 // 2) &&
           !(arguments(a)[1] isa Number)
        a = term(sqrt, arguments(a)[1])
    end

    if iscall(a) &&
       (operation(a) == ^) &&
       arguments(a)[2] isa Rational &&
       isequal((arguments(a)[2]).num, 1)
        value = demote_rational(arguments(a)[1])
        root = (arguments(a)[2]).den

        if value isa Integer && value > 0
            if isinteger(value^(1.0 / root))
                return Integer(value^(1.0 / root))
            else#find largest divisible perfect power
                outer_val = 1
                i = 2
                while i^root <= div(value, 2)
                    perfect_power = i^root
                    if value % perfect_power == 0
                        outer_val *= i
                        value = div(value, perfect_power)
                        i = 2
                        continue
                    end
                    i = i + 1
                end
                return isequal(root, 2) ? outer_val * term(sqrt, value) :
                       outer_val * term(^, value, 1 // root)
            end
        end

    end

    return a
end

#=
ex demote_rational(2//1) returns 2
ex demote_rational(2//3) returns 2//3
=#
function demote_rational(a)
    if a isa Rational && isinteger(a)
        return Integer(a)
    end
    return a
end

#=
if in quadratic form returns solutions
=#
function solve_quadratic(eq::Equation, var, single_solution)

    !iscall(eq.lhs) && return eq#make sure left side is tree form

    op = operation(eq.lhs)


    if (op == +) && isequal(degree(eq.lhs, var), 2)
        coeffs = polynomial_coeffs(eq.lhs, [var])
        a = coeffs[1][var^2]
        b = haskey(coeffs[1], var) ? coeffs[1][var] : 0
        c = coeffs[2] + coeffs[1][1] - eq.rhs

        if !(
            SymbolicUtils._occursin(var, a) ||
            SymbolicUtils._occursin(var, b) ||
            SymbolicUtils._occursin(var, c) ||
            isequal(b, 0)
        )#make sure variable in not in a b or c and that b is not zero


            sqrtPortion = reduce_root(term(sqrt, b^2 - 4 * a * c))

            if (single_solution)
                out =
                    var ~
                        SymbolicUtils.Div(sqrtPortion, 2 * a) + SymbolicUtils.Div(-b, 2 * a)
                return demote_rational(out)
            else
                out1::Any =
                    SymbolicUtils.Div(sqrtPortion, 2 * a) + SymbolicUtils.Div(-b, 2 * a)
                out2::Any =
                    -SymbolicUtils.Div(sqrtPortion, 2 * a) + SymbolicUtils.Div(-b, 2 * a)

                out1 = demote_rational(out1)
                out2 = demote_rational(out2)

                return [var ~ out1, var ~ out2]

            end
        end

    end

    return eq
end

#reverse certain functions
function inverse_funcs(eq::Equation, var)

    !iscall(eq.lhs) && return eq#make sure left side is tree form
    op = operation(eq.lhs)

    #reverse functions
    inverseOps = Dict(
        sin => asin,
        cos => acos,
        tan => atan,
        asin => sin,
        acos => cos,
        atan => tan,
        exp => log,
        log => exp,
    )

    if haskey(inverseOps, op)
        inverseOp = inverseOps[op]
        inner = arguments(eq.lhs)[1]
        eq = inner ~ term(inverseOp, eq.rhs)
    elseif (op == sqrt)
        inner = arguments(eq.lhs)[1]
        eq = inner ~ (eq.rhs)^2
    elseif (op == lambertw)
        inner = arguments(eq.lhs)[1]
        eq = inner ~ eq.rhs * term(exp, eq.rhs)
    end

    return eq
end

#solves for powers
function reverse_powers(eq::Equation, var, single_solution)
    !iscall(eq.lhs) && return eq#make sure left side is tree form
    op = operation(eq.lhs)

    if (op == ^)
        pow = eq.lhs

        pow_base = get_base(pow)
        pow_exp = get_exp(pow)

        baseHasVar = SymbolicUtils._occursin(var, pow_base)
        expoHasVar = SymbolicUtils._occursin(var, pow_exp)

        if (baseHasVar && !expoHasVar)#x^a
            twoSolutions = !single_solution && isequal(pow_exp % 2, zero(pow_exp))
            if (twoSolutions)
                eq1 = solve_single_eq(
                    pow_base ~
                        reduce_root(term(^, eq.rhs, (SymbolicUtils.Div(1, pow_exp)))),
                    var,
                )
                eq2 = solve_single_eq(
                    pow_base ~
                        -reduce_root(term(^, eq.rhs, (SymbolicUtils.Div(1, pow_exp)))),
                    var,
                )
                return [eq1, eq2]
            else
                eq =
                    pow_base ~ reduce_root(term(^, eq.rhs, (SymbolicUtils.Div(1, pow_exp))))
            end
        elseif (!baseHasVar && expoHasVar)#a^x
            eq = pow_exp ~ SymbolicUtils.Div(term(log, eq.rhs), term(log, pow_base))
        elseif (baseHasVar && expoHasVar)
            if isequal(pow_exp, pow_base)#lambert w strategy
                eq =
                    pow_exp ~ SymbolicUtils.Div(
                        term(log, eq.rhs),
                        term(lambertw, term(log, eq.rhs)),
                    )
            else#just log both sides
                eq = pow.exp * term(log, pow_base) ~ term(log, eq.rhs)
            end
        end

    end
    return eq
end

function convert_solutions_to_floats(solutions)
    out = Float64[]
    if solutions isa Equation
        check = substitute(solutions.rhs, [])
        if (check isa Number)
            push!(out, convert(Float64, check))
        else
            return nothing
        end
    elseif solutions isa Array{Equation}
        for solution in solutions
            check = substitute(solution.rhs, [])
            if (check isa Number)
                push!(out, convert(Float64, check))
            else
                return nothing
            end
        end
    end
    return out
end
