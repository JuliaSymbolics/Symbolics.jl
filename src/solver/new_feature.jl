function solve_eq_for_multivar(eq, s)
    args = arguments(unwrap(eq))
    term_tm = 1
    @info "" args
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
