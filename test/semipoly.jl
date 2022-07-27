using Symbolics
using Test

using Random

@variables x y z


@test_throws ArgumentError semipolynomial_form(x,[x],0)

d, r = semipolynomial_form(x, [x], 1)

@test d == Dict(x=>1)
@test r == 0

d, r = semipolynomial_form(x + sin(x) + 1 + y, [x], 1)

@test d == Dict(x=>1)
@test iszero(r - sin(x) - 1 - y)

d, r = semipolynomial_form(x^2+1+y, [x], 1)

@test isempty(d)
@test iszero(r - (x^2+1+y))

d, r = semipolynomial_form((x+2)^12, [x], 1)

@test d == Dict(x => 24576)
@test iszero(r + 24576x - (x+2)^12)

# 657
@test iszero(semilinear_form([x * cos(x) + x^2 * 3sin(x) + 4exp(x)], [x])[1])

pow_expr = 7^(3y + sin(y))
@test SymbolicUtils.ispow(Symbolics.unwrap(pow_expr))
dict, nl = semipolynomial_form(pow_expr, [y], Inf)
@test isempty(dict)
# expect a SymbolicUtils.Pow object instead of SymbolicUtils.Term with f = ^
@test isequal(nl, pow_expr)
@test SymbolicUtils.ispow(nl)

mul_expr = 9 * 7^y * sin(x)^4 * tan(x + y)^3
@test SymbolicUtils.ismul(Symbolics.unwrap(mul_expr))
dict, nl = semipolynomial_form(mul_expr, [x, y], Inf)
@test isempty(dict)
# expect a SymbolicUtils.Mul object instead of SymbolicUtils.Term with f = *
@test isequal(nl, mul_expr)
@test SymbolicUtils.ismul(nl)

div_expr = x / y
@test SymbolicUtils.isdiv(Symbolics.unwrap(div_expr))
dict, nl = semipolynomial_form(div_expr, [x, y], Inf)
@test isempty(dict)
# expect a SymbolicUtils.Div object instead of SymbolicUtils.Term with f = /
@test isequal(nl, div_expr)
@test SymbolicUtils.isdiv(nl)

# check negative exponent
# `SymbolicUtils.Pow` object with a negative exponent normally cannot be created.
# For example, `x^-1` returns `1 / x` which is a `SymbolicUtils.Div`.
# But `expand_derivatives` may result in a `SymbolicUtils.Pow` with a negative exponent.
sqrt_x = sqrt(x)
deriv = expand_derivatives(Differential(x)(sqrt_x)) # (1//2)*(sqrt(x)^-1)
expr = substitute(deriv, Dict(sqrt_x => y)) # (1//2)*(y^-1)
d, r = semipolynomial_form(expr, [x, y], Inf)
@test isempty(d)


@syms a b c

const components = [2, a, b, c, x, y, z, (1+x), (1+y)^2, z*y, z*x]

function verify(t, d, wrt, nl)
    try
        iszero(t - (isempty(d) ? nl : sum(k*v for (k, v) in d) + nl))
    catch err
        println("""Error verifying semi-pf result for $t
                 wrt = $wrt
                 d = $d
                 nl = $nl""")
        rethrow(err)
    end
end

seed = 0

function trial()
    global seed += 1
    Random.seed!(666+seed)
    n = rand(2:length(components)-1)
    l, r = rand(components, n), vcat(rand(components, n), [1, 0, 1//2, 1.5])

    t = *(map(1:rand(1:3)) do _
        pow = rand([1,1,1,1,1,1,1,1,2,3])
        nterms = rand(2:5)
        sum(rand(l, nterms) .* rand(r, nterms)) ^ pow
    end...)

    @show t

    for _ = 1:4
        wrt = unique(rand([a,b,c,x,y,z], rand(1:6)))
        for deg=Any[1,2,3,4,Inf]
            if deg == 1
                A, c = semilinear_form([t], wrt)
                res = iszero(A*wrt + c - [t])
                if !res
                    println("Semi-linear form is wrong: [$t]  w.r.t   $wrt ")
                    @show A c
                end
            elseif deg == 2
                A,B,v2, c = semiquadratic_form([t], wrt)
                res = iszero(A * wrt + B * v2 + c - [t])
                if !res
                    println("Semi-quadratic form is wrong: $t  w.r.t   $wrt")
                    @show A B v2 c
                end
            else
                if isfinite(deg)
                    d, nl = semipolynomial_form(t, wrt, deg)
                    @test all(x->Symbolics.pdegree(x) <= deg, keys(d))
                    for x in wrt
                        d2, enl = semipolynomial_form(expand(nl), wrt, Inf)
                        elim = all(x->Symbolics.pdegree(x)>deg, keys(d2))
                        if !elim
                            println("Imperfect elimination:")
                            @show t wrt deg nl expand(nl)
                        end
                        @test elim
                    end
                else
                    d, nl = polynomial_coeffs(t, wrt)
                end

                res = verify(t, d, wrt, nl)
                if !res
                    println("""Semi-poly form is wrong: $t  w.r.t   $wrt  deg=$deg
                            Result: $d + $nl""")
                end
            end

            @test res
        end
    end
end

for i=1:20
    @testset "fuzz semi-polynomial-form ($i/20)" begin
        trial()
    end
end
