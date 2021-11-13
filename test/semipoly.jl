using Symbolics
using Test

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


@syms a b c

const components = [2, a, b, c, x, y, z, x^2, x*y, y^2, z*y, z*x, z^2]

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

function trial()
    n = rand(2:length(components)-1)
    l, r = rand(components, n), vcat(rand(components, n), [1, 0, 1//2, 1.5])

    t = *(map(1:rand(1:3)) do _
        pow = rand([1,1,1,1,1,1,1,1,2,3])
        nterms = rand(3:20)
        sum(rand(l, nterms) .* rand(r, nterms)) ^ pow
    end...)


    for _ = 1:4
        wrt = unique(rand([a,b,c,x,y,z], rand(1:6)))
        for i=[1,2,3,4,Inf]
            if i == 1
                A, c = semilinear_form([t], wrt)
                res = iszero(A*wrt + c - [t])
                if !res
                    println("Semi-linear form is wrong: [$t]  w.r.t   $wrt ")
                    @show A c
                end
            elseif i == 2
                A,B,v2, c = semiquadratic_form([t], wrt)
                res = iszero(A * wrt + B * v2 + c - [t])
                if !res
                    println("Semi-quadratic form is wrong: $t  w.r.t   $wrt")
                    @show A B v2 c
                end
            else
                if isfinite(I)
                    d, nl = semipolynomial_form(t, wrt, i)
                else
                    d, nl = polynomial_coeffs(t, wrt)
                end

                res = verify(t, d, wrt, nl)
                if !res
                    println("""Semi-poly form is wrong: $t  w.r.t   $wrt  deg=$i
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
