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

# In the solver we might want to return Symbolics.term(sqrt, x)
# as a root. If x is a negative real number, the evaluation of such
# a root will fail since sqrt is not defined for Real n < 0.
# Thus we use Symbolics.term(Symbolics.ssqrt, x) that handles negative numbers.
# Similarly, scbrt and slog handle such similar cases.
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
        return (n)^(1 / 3)
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

const RootsOf = (SymbolicUtils.@syms roots_of(poly,var))[1]

Base.show(io::IO, f::typeof(ssqrt)) = print(io, "√")
Base.show(io::IO, r::typeof(scbrt)) = print(io, "∛")
Base.show(io::IO, r::typeof(slog)) = print(io, "slog")

function check_expr_validity(expr)
    type_expr = typeof(expr)
    valid_type = false

    if type_expr == Num || type_expr == SymbolicUtils.BasicSymbolic{Real} ||
       type_expr == Complex{Num} || type_expr == ComplexTerm{Real}
        valid_type = true
    end
    iscall(unwrap(expr)) && @assert !hasderiv(unwrap(expr)) "Differential equations are not currently supported"
    @assert valid_type "Invalid input" 
    valid_type && return nothing
    @assert isequal(expr, 0) "Invalid input"
end
function check_x(x)
    iscall(unwrap(x)) && @assert !hasderiv(unwrap(x)) "Differential equations are not currently supported"
    @assert is_singleton(unwrap(x)) "Expected a variable, got $x"
end


function check_poly_inunivar(poly, var)
    subs, filtered = filter_poly(poly, var)
    coeffs, constant = polynomial_coeffs(filtered, [var])
    return isequal(constant, 0)
end

# converts everything to BIG
function bigify(n)
    n = unwrap(n)
    if n isa ComplexTerm || n isa Float64 || n isa Irrational
        return n
    end

    if n isa SymbolicUtils.BasicSymbolic
        !iscall(n) && return n
        args = arguments(n)
        for i in eachindex(args)
            args[i] = bigify(args[i])
        end
        return n
    end

    if n isa Integer
        n = BigInt(n)
        return n
    end

    if n isa Complex
        real_part = bigify(n.re)
        im_part = bigify(n.im)
        return real_part + im_part * im
    end

    if n isa Rational && n isa Real
        n = big(n)
        return n
    end

    return n
end

function comp_rational(x, y)
    x, y = wrap(bigify(x)), wrap(bigify(y))
    if !(unwrap(x) isa AbstractFloat || x isa Complex) &&
       !(unwrap(y) isa AbstractFloat || y isa Complex)
        r = x // y
        return r
    end

    x, y = unwrap(x), unwrap(y)
    r = nothing
    if x isa ComplexF64
        real_p = real(x)
        imag_p = imag(x)
        r = Rational(real_p) // y
        if !isequal(imag_p, 0)
            r += (Rational(imag_p) // y) * im
        end
    elseif x isa Float64
        r = Rational{BigInt}(x) // y
    end

    return isequal(r, nothing) ? x / y : r
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
