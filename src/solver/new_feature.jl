function solve_interms_ofvar(eq, s; dropmultiplicity=true, warns=true)
    vars = Symbolics.get_variables(eq)
    vars = filter(v -> !isequal(v, s), vars)
    vars = wrap.(vars)

    term_tm = 1
    done = false

    # do this until no / are present
    while !done
        args = arguments(unwrap(eq))
        done = true

        for arg in args
            !iscall(arg) && continue
            oper = operation(arg)
            if oper == (/) 
                done = false
                term_tm *= wrap(arguments(arg)[2])
            end
        end
        args = [arg*term_tm for arg in args]
        eq = expand(Symbolics.term((+), unwrap.(args)...))
        term_tm = 1
    end
    coeffs, constant = polynomial_coeffs(eq, [s])
    eqs = wrap.(collect(values(coeffs)))

    symbolic_solve(eqs, vars)

end


# an attempt at using ia_solve recursively.
function find_v(eqs, v, vars)
    vars = filter(var -> !isequal(var, v), vars)
    n_eqs = deepcopy(eqs)

    if isequal(n_eqs[1], 0)
        n_eqs = n_eqs[2:end]
    end
    
    present_vars = Symbolics.get_variables(n_eqs[1])
    var = present_vars[1]
    if length(present_vars) == 1 && isequal(present_vars[1], v)
        return ia_solve(n_eqs[1], var)
    end
    if isequal(present_vars[1], v)
        var = present_vars[2]
    end

    @info "" n_eqs[1]
    @info "" var
    sols = ia_solve(n_eqs[1], var)
    isequal(sols, nothing) && return []

    for s in sols
        for j in 2:length(n_eqs)
            n_eqs[j] = substitute(n_eqs[j], Dict(var => s))
        end
        # other solutions are cut out here
        return find_v(n_eqs[2:end], v, vars)
    end

end
