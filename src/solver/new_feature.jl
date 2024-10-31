function solve_eq_for_multivar(eq, s)
    args = arguments(unwrap(eq))
    vars = Symbolics.get_variables(eq)
    vars = filter(v -> !isequal(v, s), vars)

    term_tm = 1
    if operation(unwrap(eq)) != (/)
    end

    for arg in args
        oper = operation(arg)
        if oper == (/) 
            term_tm *= wrap(arguments(arg)[2])
            @info "" term_tm
        end
    end
    args = [arg*term_tm for arg in args]
    eq = expand(Symbolics.term((+), unwrap.(args)...))
    coeffs, constant = polynomial_coeffs(eq, [s])
    eqs = wrap.(collect(values(coeffs)))

end

function solve_system_using_ia(eqs, vars, sols)
    cur_sols = []
    if isequal(eqs, [])
        return sols
    end
    try
        grb_eqs = Symbolics.groebner_basis(eqs)
        eqs = grb_eqs
    catch e
    end

    if length(new_eqs) < length(vars)
        throw("Infinite number of solutions")
    end


    eq = lastindex(eqs)
    present_vs = Symbolics.get_variables(eq)
    isequal(present_vs, []) && !isequal(eq, 0) && throw("System of equations has a contradiction.")

    if length(present_vs) == 1
        new_sols = ia_solve(eq, present_vs[1])
        append!(cur_sols, [Dict{Num, Any}(var => sol) for sol in new_sols])
        deleteat!(eqs, lastindex(eqs))

        for s in cur_sols
        end
    end
    
end
