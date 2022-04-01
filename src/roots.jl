
using PolynomialFactors
PF = PolynomialFactors
AA = PolynomialFactors.AbstractAlgebra

debug() = true

########################
##### Root Finding #####
########################

# maximal degree of polynomials that we try to solve by direct formula
const _MAX_DEGREE_DIRECT = 2

_not_implemented_sorry(ex) = throw(DomainError(ex, "Root finding not implemented, sorry."))


"""
roots(f)

Finds roots of expression `f` symbolically if possible.

"""
function roots(f)
    roots(f, get_variables(f))
end

"""
roots(f, vars)

Finds roots of expression `f` in `vars` symbolically if possible.

"""
function roots(f, vars)
    # f is a single expression here
    #=
        Currently only polynomial solving is implemented
    =#

    # transform to a convenient representation (really?)
    T = symtype(f)
    poly, p2sym, sym2p = symbol_to_poly(f)

    # construct a function to substitute polynomials back to symbols
    tosym = p -> poly_to_symbol(p, p2sym, sym2p, T)

    # !check() && throw()

    ans = Dict{Any, Any}()

    if isempty(vars)
        # pass
    elseif length(vars) == 1
        var  = only(vars)
        # what if no such var in f?
        pvar = p2sym(var)
        ans[var] = []
        _univariate_poly_roots!(poly, var, pvar, tosym, ans)
    else
        _not_implemented_sorry(f)
    end

    return ans
end

#=
    find roots of univariate polynomial
=#
function _univariate_poly_roots!(poly, var, pvar, tosym, ans)
    # var is a single variable here

    presentvars = DP.variables(poly)
    maxdeg = maximum(m -> DP.degree(m, pvar), DP.monomials(poly))

    # may want to consider this case more carefully
    # !(pvar in presentvars) && error("There is no variable $pvar in $poly")

    # if var is the only present variable
    # we can hope to output real numbers
    if length(presentvars) == 1
        _univariate_poly_roots_factor!(poly, var, pvar, maxdeg, tosym, ans)
    # otherwise solve equations directly
    # We do not want to solve degree 4..,
    # degree 3 is discussible
    elseif maxdeg <= _MAX_DEGREE_DIRECT
        _univariate_poly_roots_direct!(poly, var, pvar, maxdeg, tosym, ans)
    else
        # here is a place for non-systematic but powerful heuristics,
        # see, e.g., for cubics,
        #   http://rfdz.ph-noe.ac.at/fileadmin/Mathematik_Uploads/ACDCA/TIME2016/Beaudin_ao_Third_Degree_Polynomial_Equations_-_Presentation.pdf
        #   https://www.shsu.edu/kws006/professional/Concepts_files/SolvingCubics.pdf
        _not_implemented_sorry(poly)
    end
end

# hack,
# return symbolic coefficient of term t in monomial x,
function extract_coefficient(t, x)
    monom = DP.leadingmonomial(DP.div(t, x))
    coeff = DP.coefficient(t)
    @info " " coeff*monom
    coeff * monom
end

# return a, b assuming poly is ax - b
function decompose_linear(poly, pvar)
    ax = argmax(m -> DP.degree(m, pvar), DP.terms(poly))
    a = extract_coefficient(ax, pvar^DP.degree(ax, pvar))
    b  = -poly + ax
    a, b
end

# return a, b, c assuming poly is ax^2 + bx + c
function decompose_quadratic(poly, pvar)
    ax = argmax(m -> DP.degree(m, pvar), DP.terms(poly))
    a = extract_coefficient(ax, pvar^DP.degree(ax, pvar))
    poly = poly - ax
    bx = argmax(m -> DP.degree(m, pvar), DP.terms(poly))
    if DP.degree(bx, pvar) == 0
        b = 0
        c = poly
    else
        b = extract_coefficient(bx, pvar^DP.degree(bx, pvar))
        c = poly - bx
    end
    a, b, c
end

# this is so unstable yay
_simplesqrt(x)          = term(sqrt, x)
_simplesqrt(x::Complex) = term(sqrt, x)
function _simplesqrt(x::Number)
    if floor(sqrt(x)^2) == x
        floor(typeof(x), sqrt(x))
    else
        term(sqrt, x)
    end
end

#=
    find roots of univariate polynomial of small degree directly
=#
function _univariate_poly_roots_direct!(poly, var, pvar, maxdeg, tosym, ans)
    if debug()
        @info "direct"
    end

    if maxdeg == 0
        # pass
    # if linear
    elseif maxdeg == 1
        # assuming (not unreasonable)
        # all terms containing var were condensed into one
        a, b = decompose_linear(poly, pvar)
        root = tosym(b) // tosym(a)
        ans[var] = root
    # if quadratic
    elseif maxdeg == 2
        a, b, c = decompose_quadratic(poly, pvar)
        d = b^2 - 4a*c
        if DP.isconstant(d) && sign(DP.leadingcoefficient(d)) < 0
            d = Complex(DP.leadingcoefficient(d))
        end
        c, b, a, d = tosym(c), tosym(b), tosym(a), tosym(d)
        if debug()
            println("c = $c, b = $b, a = $a, d = $d")
        end
        d = _simplesqrt(d)
        root1 = (-b - d) // (2a)
        root2 = (-b + d) // (2a)
        push!(ans[var], simplify(root1))
        push!(ans[var], simplify(root2))
    else
        _not_implemented_sorry(poly)
    end
end

#=
    find roots of univariate polynomial by factorization
=#
function _univariate_poly_roots_factor!(poly, var, pvar, maxdeg, tosym, ans)
    if debug()
        @info "factor"
    end

    tobefactorized = PF.as_poly(dp_densecoeffs(poly, pvar, maxdeg))
    factors = PF.poly_factor(tobefactorized)
    # drop multiplicities here for now
    debug() && println(factors)
    for (factor, mult) in factors
        deg = AA.degree(factor)
        # if trivial
        if deg == 0
            # pass
        # if linear
        elseif deg == 1
            root = -AA.constant_coefficient(factor)
            push!(ans[var], root)
        elseif deg <= _MAX_DEGREE_DIRECT
            dpfactor = aa_to_dp(factor, pvar)
            _univariate_poly_roots_direct!(dpfactor, var, pvar, deg, tosym, ans)
        else
            _not_implemented_sorry(poly)
        end
    end
end

#  =(
#   sad
dp_densecoeffs(f, x, deg) = [DP.coefficient(f, x^i) for i in 0:deg]
aa_densecoeffs(f, x, deg) = [AA.coeff(f, i) for i in 0:deg]

# abstract algebra poly to dynamic polynomial
function aa_to_dp(f, x)
    fcoeffs = aa_densecoeffs(f, AA.var(AA.parent(f)), AA.degree(f))
    mapreduce(ic -> ic[2]*x^(ic[1]-1), +, enumerate(fcoeffs))
end

########################
#### System Solving ####
########################

"""
roots(system)

Finds roots of system `system` symbolically if possible.

"""
function roots(system::AbstractArray)
    roots(system, get_variables(system))
end

function roots(system::AbstractArray, vars)
    ans = Dict{Any, Any}()

    # transform to a convenient representation (?)
    T = symtype(first(system))
    polys, p2sym, sym2p = symbol_to_poly(system)

    tosym = p -> poly_to_symbol(p, p2sym, sym2p, T)

    # the sort would not be needed
    basis = sort(Groebner.groebner(polys), by=DP.leadingmonomial, rev=true)
    pvars = sort(DP.variables(first(basis)), rev=false)

    npvars = length(pvars)

    if debug()
        println(basis)
        println(pvars)
    end

    if npvars == 1
        # _multivatiate1_system_roots!(basis, var, pvar, tosym, ans)
    elseif npvars == 2
        _multivatiate2_system_roots!(basis, vars, pvars, tosym, ans)
    elseif npvars == 3
        # _multivatiate3_system_roots!(basis, vars, pvars, tosym, ans)
    else
        _not_implemented_sorry(system)
    end

    ans
end

function _multivatiate2_system_roots!(basis, vars, pvars, tosym, ans)

    eq1   = basis[1]
    pvar1 = pvars[1]
    var1  = tosym(pvar1)
    deg1   = DP.degree(DP.leadingmonomial(eq1), pvar1)

    anslocal1 = Dict{Any, Any}()
    anslocal1[var1] = []
    _univariate_poly_roots_direct!(eq1, var1, pvar1, deg1, tosym, anslocal1)

    if debug()
        println("local1: $anslocal1")
    end

    for root1 in anslocal1[var1]
        specbasis = specializesystem(basis, pvar1, root1)
        if debug()
            println("spec: $specbasis")
        end

        eq2   = specbasis[2]
        pvar2 = pvars[2]
        var2  = tosym(pvar2)
        deg2   = DP.degree(DP.leadingmonomial(eq2), pvar2)

        if !haskey(ans, (var1, var2))
            ans[(var1, var2)] = []
        end

        anslocal2 = Dict{Any, Any}()
        anslocal2[var2] = []
        _univariate_poly_roots_direct!(eq2, var2, pvar2, deg2, tosym, anslocal2)

        for root2 in anslocal2[var2]
            push!(ans[(var1, var2)], (root1, root2))
        end
    end
end

# substitute variable `var` with `root` value everywhere in `system`
function specializesystem(system, var, root)
    specsystem = similar(system)
    for (i, f) in enumerate(system)
        specsystem[i] = DP.subs(f, var=>root)
    end
    specsystem
end
