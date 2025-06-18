function cross_multiply(eq)
    og_oper = operation(unwrap(eq))
    done = true
    loop_add = false

    term_tm = 1
    if og_oper === (/)
        done = false
        args = arguments(unwrap(eq))
        eq = wrap(args[1])
    end

    # do this until no / are present
    if og_oper === (+)
        while !loop_add
            args = arguments(unwrap(eq))
            loop_add = true

            for arg in args
                !iscall(arg) && continue
                oper = operation(arg)
                if oper == (/) 
                    done = false
                    loop_add = false
                    term_tm *= wrap(arguments(arg)[2])
                end
            end
            args = [arg*term_tm for arg in args]
            eq = expand(Symbolics.term((+), unwrap.(args)...))
            term_tm = 1
        end
    end

    if done
        return eq
    else
        return cross_multiply(eq)
    end
end
"""
    solve_interms_ofvar(eq, s; dropmultiplicity=true, warns=true)
This special case solver expects a single equation in multiple variables and a
variable `s` (this can be any Num, `s` is used for convenience). The function generates
a system of equations to by observing the coefficients of the powers of `s` present in `eq`.
E.g. a system would look like `a+b = 1`, `a-2b = 3` for the eq `(a+b)s + (a-2b)s^2 - (1)s - (3)s^2 = 0`.
After generating this system, it calls `symbolic_solve`, which uses `solve_multivar`. `symbolic_solve` was chosen
instead of `solve_multivar` because it postprocesses the roots in order to simplify them and make them more user friendly.

Generation of system uses cross multiplication in order to simplify the equation and convert it
to a polynomial like shape.


# Arguments
- eq: Single symbolics Num or SymbolicUtils.BasicSymbolic. This is equated to 0 and then solved. E.g. `expr = x+2`, we solve `x+2 = 0`

- s: Variable to "isolate", i.e. ignore and generate the system of equations based on this variable's coefficients.

- dropmultiplicity (optional): Print repeated roots or not?

- warns (optional, this is not used currently): Warn user when something is wrong or not.

# Examples
```jldoctest
julia> @variables a b x s;

julia> eq = (a*x^2+b)*s^2 - 2s^2 + 2*b*s - 3*s + 2(x^2)*(s^3) + 10*s^3;

julia> Symbolics.solve_interms_ofvar(eq, s)
2-element Vector{Any}:
 Dict{Num, Any}(a => -1//10, b => 3//2, x => (0 - 1im)*√(5))
 Dict{Num, Any}(a => -1//10, b => 3//2, x => (0 + 1im)*√(5))
```
```jldoctest
julia> eq = ((s^2 + 1)/(s^2 + 2*s + 1)) - ((s^2 + a)/(b*c*s^2 + (b+c)*s + d));

julia> Symbolics.solve_interms_ofvar(eq, s)
1-element Vector{Any}:
 Dict{Num, Any}(a => 1, d => 1, b => 1, c => 1)
```
"""
function solve_interms_ofvar(eq, s; dropmultiplicity=true, warns=true)
    @assert iscall(unwrap(eq))
    vars = Symbolics.get_variables(eq)
    vars = filter(v -> !isequal(v, s), vars)
    vars = wrap.(vars)

    term_tm = 1
        
    eq = cross_multiply(eq)
    coeffs, constant = polynomial_coeffs(eq, [s])
    eqs = wrap.(collect(values(coeffs)))

    symbolic_solve(eqs, vars, dropmultiplicity=dropmultiplicity, warns=warns)
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
