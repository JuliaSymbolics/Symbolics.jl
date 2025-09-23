
"""
    sub(subs, place_to_sub)
Helper function for filter_poly. Generates a symbolics variable and adds it
to subs (changes the original array passed to the sub function). returns the
subbed value so the filter_poly function could change it in its scope.

# Arguments
- subs: Vector of dicts which consist of var_subbed => value_subbed in filter_poly
- place_to_sub: The place which should be substituted by a new variable.

"""
function sub(subs, place_to_sub)
    sub_var = gensym()
    sub_var = (@variables $sub_var)[1]

    subs[sub_var] = place_to_sub
    place_to_sub = sub_var.val

    return place_to_sub
end

"""
    clean_f(filtered_expr, var, subs)

Helper function for `filter_poly` which is called directly before returning
the `filtered_expression`. This function aims to get the filtered expressions
resulting from `filter_poly` ready to get used by solve, `get_roots`, and any other
function in the library. An important feature of it is that it simplifies fractions
to make the output "valid" in the eyes of `factor_use_nemo` and other Nemo functions.

# Arguments
- `filtered_expr`: The output from `_filter_poly`.
- var: The variable which is filtered for.
- subs: Vector of dicts which consist of `var_subbed` => `value_subbed` in `_filter_poly`.

"""
function clean_f(filtered_expr, var, subs)
    filtered_expr = simplify_fractions(simplify(real(filtered_expr)))
    unwrapped_f = unwrap(filtered_expr)
    !iscall(unwrapped_f) && return filtered_expr
    oper = operation(unwrapped_f)
    assumptions = []

    if oper === (/)
        if !all(isequal(var, x) for x in get_variables(denominator(unwrapped_f)))
            filtered_expr = numerator(unwrapped_f)
            push!(assumptions, substitute(denominator(unwrapped_f), subs, fold=false))
        end
    end
    return filtered_expr, assumptions
end

"""
    filter_stuff(expr)
Helper function for _filter_poly and filter_poly which aims to filter
constants and things that dont contain variables like term(sqrt, 2).
Filters only if the passed expr is scary (i.e. not a Rational or Int).

# Arguments
- expr: The detected constant in _filter_poly

# Examples
```jldoctest
julia> filter_stuff(term(sqrt, 2))
(Dict{Any, Any}(var"##278" => sqrt(2)), var"##278")

julia> filter_stuff(123)
(Dict{Any, Any}(), 123)
```
"""
function filter_stuff(expr)
    expr = value(expr)
    if expr isa Integer
        return Dict(), expr
    elseif expr isa Rational || expr isa AbstractFloat || expr isa Complex
        return Dict(), comp_rational(expr, 1)
    else
        expr = isequal(expr, true) ? 1 : expr
        subs = Dict{Any, Any}()

        expr = sub(subs, expr)
        return subs, expr
    end
end

"""
    _filter_poly(expr, var)

Main mechanism for filter_poly. Traverses arguments as deep 
as needed/makes sense by recalling itself as many times as necessary.
Output is then returned to filter_poly. This function should not be used
(although functional) instead of filter_poly since it does not contain
the final touch of clean_f, which is needed if the output is going to get passed
to other functions.

# Arguments
- expr: Expression to be filtered passed from filter_poly
- var: Var that the expression should be filtered with respect to.

# Examples
```jldoctest
julia> _filter_poly(x + sqrt(2), x)
(Dict{Any, Any}(var"##239" => 1.4142135623730951), var"##239" + x)

julia> RootFinding._filter_poly(x*sqrt(2), x)
(Dict{Any, Any}(var"##240" => 1.4142135623730951), var"##240"*x)
```
"""
function _filter_poly(expr, var)
    expr = unwrap(expr)
    vars = get_variables(expr)
    if !isempty(vars) && isequal(first(vars), expr)
        return (Dict{Any, Any}(), expr)
    elseif isempty(vars)
        return filter_stuff(expr)
    end

    args = copy(parent(arguments(expr)))
    if symtype(expr) <: Complex
        subs1, subs2 = Dict(), Dict()
        expr1, expr2 = 0, 0
        rr = real(expr)
        ii = imag(expr)
        if !isequal(rr, 0)
            subs1, expr1 = _filter_poly(rr, var)
        end
        if !isequal(ii, 0)
            subs2, expr2 = _filter_poly(ii, var)
        end

        subs = merge(subs1, subs2)
        i_var = gensym()
        i_var = (@variables $i_var)[1]

        subs[i_var] = im
        expr = unwrap(expr1 + i_var * expr2)

        args = map(unwrap, arguments(expr))
        oper = operation(expr)
        return subs, term(oper, args...)
    end
    subs = Dict{Any, Any}()
    for (i, arg) in enumerate(args)
        # handle constants
        arg = value(arg)
        vars = get_variables(arg)
        if isempty(vars)
            if arg isa Integer
                args[i] = Const{VartypeT}(bigify(args[i]))
                continue
            elseif arg isa Rational || arg isa AbstractFloat || arg isa Complex
                args[i] = Const{VartypeT}(comp_rational(arg, 1))
                continue
            end
            args[i] = Const{VartypeT}(sub(subs, args[i]))
            continue
        end

        # handle "x" as an argument
        if length(vars) == 1
            if isequal(arg, var) || isequal(first(vars), arg)
                continue
            end
        end

        oper = operation(arg)
        monomial = copy(parent(arguments(arg)))
        if oper === (^)
            if any(arg -> isequal(arg, var), monomial)
                continue
            end
            # filter(args[1]), filter[args[2]] and then merge
            subs1, __monomial_1 = _filter_poly(monomial[1], var)
            subs2, __monomial_2 = _filter_poly(monomial[2], var)
            monomial[1] = Const{VartypeT}(__monomial_1)
            monomial[2] = Const{VartypeT}(__monomial_2)

            merge!(subs, merge(subs1, subs2))
            args[i] = maketerm(typeof(arg), oper, monomial, metadata(arg))
            continue
        end

        if oper === (*)
            subs_of_monom = Dict{Any, Any}()
            for (j, x) in enumerate(monomial)
                vars = get_variables(x)
                if (!isempty(vars) && isequal(first(vars), x))
                    continue
                elseif x isa Integer
                    monomial[j] = bigify(monomial[j])
                    continue
                elseif x isa Rational || x isa AbstractFloat || x isa Complex
                    monomial[j] = comp_rational(x, 1)
                    continue
                end
                # filter each arg and then merge
                new_subs, __monomial_j = _filter_poly(monomial[j], var)
                monomial[j] = Const{VartypeT}(__monomial_j)
                merge!(subs_of_monom, new_subs)
            end
            merge!(subs, subs_of_monom)
            args[i] = maketerm(typeof(arg), oper, monomial, metadata(arg))
            continue
        end

        if oper === (/) || oper === (+)
            for (j, x) in enumerate(monomial)
                new_subs, new_filtered = _filter_poly(monomial[j], var)
                merge!(subs, new_subs)
            end
            continue
        end
    end

    args = map(unwrap, args)
    oper = operation(expr)
    expr = maketerm(typeof(expr), oper, args, metadata(expr))
    return subs, expr
end

"""
    filter_poly(og_expr, var)

Filters polynomial expressions from weird irrationals that mess up
the functionality of other helper functions for our solvers such as
polynomial_coeffs and factor_use_nemo.

# Arguments
- og_expr: Original passed expression. This is deepcopied so that the
user's expression is not altered unwillingly.
- var: Var that the expression should be filtered with respect to.

# Examples
```jldoctest
julia> filter_poly(x + 2im, x)
(Dict{Any, Any}(var"##244" => im), 2var"##244" + x)

julia> filter_poly((1/im)*x + 3*y*z, x)
(Dict{Any, Any}(var"##245" => -1.0, var"##246" => im), 3y*z + var"##245"*var"##246"*x)

julia> filter_poly((x+1)*term(log, 3), x)
(Dict{Any, Any}(var"##247" => log(3)), var"##247"*(1 + x))
```
"""
function filter_poly(og_expr, var; assumptions=false)
    expr = deepcopy(og_expr)
    expr = unwrap(expr)
    vars = get_variables(expr)

    # handle edge cases
    if !isempty(vars) && isequal(first(vars), expr)
        assumptions && return Dict{Any, Any}(), expr, []
        return (Dict{Any, Any}(), expr)
    elseif isempty(vars)
        assumptions && return filter_stuff(expr), []
        return filter_stuff(expr)
    end

    # core filter
    subs, expr = _filter_poly(expr, var)

    # reassemble expr to avoid variables remembering original values issue and clean
    args = arguments(expr)
    oper = operation(expr)
    new_expr, assum_array = clean_f(term(oper, args...; type = symtype(expr), shape = shape(expr)), var, subs)

    assumptions && return subs, new_expr, assum_array
    return subs, new_expr
end

function filter_poly(og_expr; assumptions=false)
    new_var = gensym()
    new_var = (@variables $(new_var))[1]
    return filter_poly(og_expr, new_var; assumptions=assumptions)
end


"""
    sdegree(coeffs, var)
Gets the degree of a polynomial by traversing the `coeffs`
output from polynomial_coeffs.

# Arguments
- coeffs: output from polynomial_coeffs
- var: var present in polynomial_coeffs

# Examples
```jldoctest
julia> coeffs, constant = polynomial_coeffs(x^2 + x + 1, [x])
(Dict{Any, Any}(x^2 => 1, x => 1, 1 => 1), 0)

julia> sdegree(coeffs, x)
2


julia> coeffs, constant = polynomial_coeffs(x^12 + 3x^2, [x])
(Dict{Any, Any}(x^2 => 3, x^12 => 1), 0)

julia> sdegree(coeffs, x)
12
```
"""
function sdegree(coeffs, var)
    degree = 0
    vars = collect(keys(coeffs))
    for n in vars
        SymbolicUtils._isone(n) && continue
        isequal(n, var) && degree > 1 && continue

        if isequal(n, var) && degree < 1
            degree = 1
            continue
        end

        args = arguments(n)
        degree = max(unwrap_const(args[2]), degree)
    end
    return degree
end
