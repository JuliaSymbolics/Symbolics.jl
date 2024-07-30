# general helpers
function ssubs(expr, dict)
    if haskey(dict, expr)
        return dict[expr]
    end
    expr = unwrap(expr)
    if iscall(expr)
        op = ssubs(operation(expr), dict)
        args = map(arguments(expr)) do x
            ssubs(x, dict)
        end
        op(wrap.(args)...)
    else
        return expr
    end
end

function ssqrt(n)
    n = unwrap(n)

    if n isa Real
        n > 0 && return sqrt(n)
        return sqrt(complex(n))
    end

    if n isa Complex 
        return sqrt(n)
    end

    if n isa SymbolicUtils.BasicSymbolic{Real}
        return term(ssqrt, n)
    end

end


function scbrt(n)
    n = unwrap(n)

    if n isa Real
        return cbrt(n)
    end

    if n isa Complex
        return (n)^(1/3)
    end

    if n isa SymbolicUtils.BasicSymbolic{Real}
        return term(scbrt, n)
    end

end

function slog(n)
    n = unwrap(n)

    if n isa Real
        return n > 0 ? log(n) : log(complex(n))
    end

    if n isa Complex
        return log(n)
    end

    return term(slog, n)
end


struct RootsOf
    poly::Num
    var::Num
end

Base.show(io::IO, r::RootsOf) = print(io, "roots_of(", r.poly, ", ", x, ")")
Base.show(io::IO, f::typeof(ssqrt)) = print(io, "√")
Base.show(io::IO, r::typeof(scbrt)) = print(io, "∛")

# not sure if this is a good idea as it can hide from the
# user when it misbehaves
# Base.show(io::IO, r::typeof(slog)) = print(io, "log")


function check_expr_validity(expr)
    type_expr = typeof(expr)
    valid_type = false

    if type_expr == Num || type_expr == SymbolicUtils.BasicSymbolic{Real} || 
    type_expr == Complex{Num} || type_expr == ComplexTerm{Real}
        valid_type = true
    end
    @assert valid_type && return true
    @assert isequal(expr, 0) "Invalid input"
end

function check_poly_inunivar(poly, var)
    subs, filtered = filter_poly(poly, var)
    coeffs, constant = polynomial_coeffs(filtered, [var])
    return isequal(constant, 0)
end

function f_numbers(n)
    n = unwrap(n)
    if n isa ComplexTerm || n isa Float64 || n isa Irrational
        return n
    end

    if n isa SymbolicUtils.BasicSymbolic
        !iscall(n) && return n
        args = arguments(n)
        for i in eachindex(args)
            args[i] = f_numbers(args[i])
        end
        return n
    end

    if n isa Integer
        n = BigInt(n)
        return n
    end

    if n isa Complex
        real_part = f_numbers(n.re)
        im_part = f_numbers(n.im)
        return real_part + im_part*im
    end

    if n isa Rational && n isa Real
        n = big(n)
        return n
    end

    return n
end

function comp_rational(x,y)
    x, y  = f_numbers(x), f_numbers(y)
    try
        r = x//y
        return r
    catch e
        return x/y
    end
end



### multivar stuff ###
function contains_var(var, vars)
    for variable in vars
        if isequal(var, variable)
            return true
        end
    end
    return false
end

function add_sol!(solutions, new_sols, var, index)
    sol_used = solutions[index]
    deleteat!(solutions, index)
    for new_sol in new_sols
        sol_used[var] = new_sol
        push!(solutions, deepcopy(sol_used))
    end
    return solutions
end

function add_sol_to_all(solutions, new_sols, var)
    existing_solutions = deepcopy(solutions)
    solutions = []
    for new_sol in new_sols
        copy_sol = deepcopy(existing_solutions)
        for i in eachindex(copy_sol)
            copy_sol[i][var] = new_sol
        end
        append!(solutions, copy_sol)
    end
    return solutions
end

