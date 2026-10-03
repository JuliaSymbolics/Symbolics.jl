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
        isnan(n) && return n
        n > 0 && return sqrt(n)
        return sqrt(complex(n))
    end

    if n isa Complex
        return sqrt(n)
    end

    if symtype(n) === Real
        return term(ssqrt, n)
    end
end

SymbolicUtils.promote_symtype(::typeof(ssqrt), ::Type{T}) where {T} = T
SymbolicUtils.promote_shape(::typeof(ssqrt), @nospecialize(sh::SymbolicUtils.ShapeT)) = sh

@register_derivative ssqrt(x) I begin
    substitute(@derivative_rule(sqrt(x), I), sqrt => ssqrt)
end

function scbrt(n)
    n = unwrap(n)

    if n isa Real
    isnan(n) && return n
        return cbrt(n)
    end

    if n isa Complex
        return (n)^(1 / 3)
    end

    if symtype(n) === Real
        return term(scbrt, n)
    end
end

SymbolicUtils.promote_symtype(::typeof(scbrt), ::Type{T}) where {T} = T
SymbolicUtils.promote_shape(::typeof(scbrt), @nospecialize(sh::SymbolicUtils.ShapeT)) = sh
@register_derivative scbrt(x) I begin
    substitute(@derivative_rule(cbrt(x), I), cbrt => scbrt)
end

function slog(n)
    n = unwrap(n)

    if n isa Real
        isnan(n) && return n
        return n > 0 ? log(n) : log(complex(n))
    end

    if n isa Complex
        return log(n)
    end

    return term(slog, n)
end

SymbolicUtils.promote_symtype(::typeof(slog), ::Type{T}) where {T} = T
SymbolicUtils.promote_shape(::typeof(slog), @nospecialize(sh::SymbolicUtils.ShapeT)) = sh

@register_derivative slog(x) I begin
    substitute(@derivative_rule(log(x), I), log => slog)
end

const RootsOf = (SymbolicUtils.@syms roots_of(poly,var))[1]

Base.show(io::IO, f::typeof(ssqrt)) = print(io, "√")
Base.show(io::IO, r::typeof(scbrt)) = print(io, "∛")
Base.show(io::IO, r::typeof(slog)) = print(io, "slog")

function check_expr_validity(expr)
    type_expr = typeof(expr)
    valid_type = false
    st = symtype(expr)
    if type_expr <: Number || type_expr == Num || st <: Real ||
       type_expr == Complex{Num} || st <: Complex{Real}
        valid_type = true
    end
    iscall(unwrap(expr)) && @assert !hasderiv(unwrap(expr)) "Differential equations are not currently supported"
    @assert valid_type "Invalid input of type $type_expr (symtype $st)"
    return valid_type && return nothing
end
function check_x(x)
    iscall(unwrap(x)) && @assert !hasderiv(unwrap(x)) "Differential equations are not currently supported"
    @assert is_singleton(unwrap(x)) "Expected a variable, got $x"
end


function check_poly_inunivar(poly, var)
    subs, filtered = filter_poly(poly, var)
    coeffs, constant = polynomial_coeffs(filtered, var isa Array ? var : [var])
    return SymbolicUtils._iszero(constant)
end

# converts everything to BIG
"""
    exact_div(x, y)

Exact division that also accepts symbolic operands.

`//` has no method for `BasicSymbolic`, so a parametric coefficient makes the
closed-form root formulas in `univar.jl` throw a `MethodError`. Dividing two
symbolic expressions is exact already, since it builds a `Div`, so this routes
through `/` unless both operands are integers or rationals — where `/` would
return a float and lose the exactness those formulas rely on.
"""
exact_div(x, y) = x / y
exact_div(x::Union{Integer, Rational}, y::Union{Integer, Rational}) = x // y

function bigify(n)
    n = value(n)
    if n isa Float64 || n isa Irrational
        return n
    end

    if n isa SymbolicUtils.BasicSymbolic
        !iscall(n) && return n
        args = copy(parent(arguments(n)))
        for i in eachindex(args)
            args[i] = Const{VartypeT}(bigify(args[i]))
        end
        n = maketerm(typeof(n), operation(n), args, metadata(n))
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
    x, y = bigify(unwrap(x)), bigify(unwrap(y))
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

function is_exact_polynomial(expr)
    expr = unwrap(expr)
    if SymbolicUtils.isconst(expr) || expr isa Number
        c = value(expr)
        return c isa Integer || (c isa Rational && !iszero(denominator(c)))
    end
    is_singleton(expr) && return true
    iscall(expr) || return false
    op = operation(expr)
    args = arguments(expr)
    if op === (+) || op === (*)
        return all(is_exact_polynomial, args)
    elseif op === (^)
        exponent = value(args[2])
        return exponent isa Integer && exponent >= 0 && is_exact_polynomial(args[1])
    end
    return false
end

canonical_affine_coefficient(expr) = simplify_fractions(expand(expr))

function exact_affine_solve(eqs, vars)
    length(eqs) == length(vars) || return nothing
    isempty(vars) && return nothing
    A, bvec, islinear = linear_expansion(wrap.(bigify.(eqs)), vars)
    islinear || return nothing
    x_set = Set(unwrap(v) for v in vars)
    for e in Iterators.flatten((A, bvec))
        any(v -> v in x_set, get_variables(e)) && return nothing
    end
    A = canonical_affine_coefficient.(wrap.(bigify.(A)))
    bvec = -wrap.(bigify.(bvec))
    n = length(vars)
    for k in 1:n
        pivot = findfirst(i -> !_iszero(A[i, k]), k:n)
        isnothing(pivot) && return nothing
        pivot += k - 1
        if pivot != k
            A[k, :], A[pivot, :] = A[pivot, :], A[k, :]
            bvec[k], bvec[pivot] = bvec[pivot], bvec[k]
        end
        for i in (k + 1):n
            _iszero(A[i, k]) && continue
            factor = canonical_affine_coefficient(A[i, k] / A[k, k])
            for j in (k + 1):n
                A[i, j] = canonical_affine_coefficient(A[i, j] - factor * A[k, j])
            end
            bvec[i] = canonical_affine_coefficient(bvec[i] - factor * bvec[k])
            A[i, k] = 0
        end
    end
    return symsub!(UpperTriangular(A), bvec)
end

# Strip outer integer powers and nonzero constant factors so that f^n and c*f^n
# share the zero set of f when multiplicities are discarded.
function drop_outer_multiplicities(expression)
    expression = unwrap(expression)
    changed = true
    while changed && iscall(expression)
        changed = false
        op = operation(expression)
        args = arguments(expression)
        if isequal(op, ^) && SymbolicUtils.isconst(args[2])
            a2 = unwrap_const(args[2])
            if a2 isa Integer && a2 > 0
                expression = unwrap(args[1])
                changed = true
                continue
            end
        elseif isequal(op, *)
            new_factors = Any[]
            local_changed = false
            for a in args
                a = unwrap(a)
                if iscall(a) && isequal(operation(a), ^)
                    aa = arguments(a)
                    if SymbolicUtils.isconst(aa[2])
                        exp = unwrap_const(aa[2])
                        if exp isa Integer && exp > 0
                            push!(new_factors, aa[1])
                            local_changed = true
                            continue
                        end
                    end
                    push!(new_factors, a)
                elseif SymbolicUtils.isconst(a) || a isa Number
                    if isequal(a, 0) || (a isa Number && iszero(a))
                        return wrap(0)
                    end
                    local_changed = true
                else
                    push!(new_factors, a)
                end
            end
            isempty(new_factors) && return wrap(1)
            expression = length(new_factors) == 1 ? unwrap(new_factors[1]) :
                unwrap(*(wrap.(new_factors)...))
            changed = local_changed
            continue
        end
    end
    return wrap(expression)
end
