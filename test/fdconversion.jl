
import Symbolics
import FastDifferentiation as FD
using FDConversion
using Random
using Test

#Function used in tests for FDConversion
function P(::Type{T}, l, m, z::T) where {T}
    if l == 0 && m == 0
        return T(1)
    elseif l == m
        return (1 - 2m) * P(T, m - 1, m - 1, z)
    elseif l == m + 1
        return (2m + 1) * z * P(T, m, m, z)
    else
        return ((2l - 1) / (l - m) * z * P(T, l - 1, m, z) - (l + m - 1) / (l - m) * P(T, l - 2, m, z))
    end
end


function S(::Type{T}, m, x::T, y::T) where {T}
    if m == 0
        return T(0)
    else
        return x * C(T, m - 1, x, y) - y * S(T, m - 1, x, y)
    end
end


function C(::Type{T}, m, x::T, y::T) where {T}
    if m == 0
        return T(1)
    else
        return x * S(T, m - 1, x, y) + y * C(T, m - 1, x, y)
    end
end

function factorial_approximation(::Type{T}, x) where {T}
    local n1 = x
    sqrt(2 * T(π) * n1) * (n1 / T(ℯ) * sqrt(n1 * sinh(1 / T(n1)) + 1 / (810 * T(n1)^6)))^n1
end


function compare_factorial_approximation()
    for n in 1:30
        println("n $n relative error $((factorial(big(n))-factorial_approximation(BigFloat,n))/factorial(big(n)))")
    end
end


function N(::Type{T}, l, m) where {T}
    @assert m >= 0
    if m == 0
        return sqrt((2l + 1 / (4 * T(π))))
    else
        # return sqrt((2l+1)/2π * factorial(big(l-m))/factorial(big(l+m)))
        #use factorial_approximation instead of factorial because the latter does not use Stirlings approximation for large n. Get error for n > 2 unless using BigInt but if use BigInt get lots of rational numbers in symbolic result.
        return sqrt((2l + 1) / 2 * T(π) * factorial_approximation(T, l - m) / factorial_approximation(T, l + m))
    end
end


"""l is the order of the spherical harmonic"""
function Y(::Type{T}, l, m, x::T, y::T, z::T) where {T}
    @assert l >= 0
    @assert abs(m) <= l
    if m < 0
        return N(T, l, abs(m)) * P(T, l, abs(m), z) * S(T, abs(m), x, y)
    else
        return N(T, l, m) * P(T, l, m, z) * C(T, m, x, y)
    end
end

Y(l, m, x::T, y::T, z::T) where {T<:FD.Node} = Y(FD.Node, l, m, x, y, z)
Y(l, m, x::T, y::T, z::T) where {T<:Number} = Y(T, l, m, x, y, z)


function SHFunctions(max_l, x, y, z)
    @assert typeof(x) == typeof(y) == typeof(z)
    result = Vector(undef, 0)

    for l in 0:max_l-1
        for m in -l:l
            push!(result, Y(l, m, x, y, z))
        end
    end

    return result
end



@testset "conversion from Symbolics to FD" begin
    Symbolics.@variables x y

    symbolics_expr = x^2 + y * (x^2)
    dag, tmp = to_fd(symbolics_expr)
    vars = collect(values(tmp))
    fdx, fdy = FD.value(vars[1]) == :x ? (vars[1], vars[2]) : (vars[2], vars[1]) #need to find the variables since they can be in any order

    correct_fun = FD.make_function([dag], [fdx, fdy])


    #verify that all the node expressions exist in the dag. Can't rely on them being in a particular order because Symbolics can
    #arbitrarily choose how to reorder trees.
    num_tests = 100
    rng = Random.Xoshiro(8392)
    for _ in 1:num_tests
        (xval, yval) = rand(rng, 2)
        FDval = correct_fun([xval, yval])[1]
        Syval = Symbolics.substitute(symbolics_expr, Dict([(x, xval), (y, yval)]))

        @test isapprox(FDval, Syval.val)

    end
end


@testset "conversion from FD to Symbolics" begin
    order = 8
    FD.@variables x y z

    FD_funcs = FD.Node.(SHFunctions(order, x, y, z))
    Sym_funcs, variables = to_symbolics(FD_funcs)

    sx, sy, sz = map(p -> variables[p], [x, y, z])

    FD_eval = FD.make_function(FD_funcs, [x, y, z])
    rng = Random.Xoshiro(8392)
    num_tests = 100
    for _ in 1:num_tests
        tx, ty, tz = rand(rng, BigFloat, 3)
        subs = Dict([sx => tx, sy => ty, sz => tz])
        res = Symbolics.substitute.(Sym_funcs, Ref(subs))

        FD_res = FD_eval([tx, ty, tz])

        @test isapprox(FD_res, res, atol=1e-12)
    end
end


