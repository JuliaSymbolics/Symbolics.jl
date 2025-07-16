using ForwardDiff
using Symbolics
using Symbolics.SymbolicUtils
using Symbolics.SymbolicUtils.SpecialFunctions
using Symbolics.NaNMath
using Test

SF = SymbolicUtils.SpecialFunctions

@variables x

# Test functions from Symbolics #
#-------------------------------#

for f ∈ SymbolicUtils.basic_monadic
    fun = eval(:(ξ ->($f)(ξ)))

    fd = ForwardDiff.derivative(fun, x)
    sym = Symbolics.Differential(x)(fun(x)) |> expand_derivatives

    @test isequal(fd, sym)
end

for f ∈ SymbolicUtils.monadic
    # The polygamma and trigamma functions seem to be missing rules in ForwardDiff.
    # The abs rule uses conditionals and cannot be used with Symbolics.Num.
    # acsc, asech, NanMath.log2 and NaNMath.log10 are tested separately
    if f ∈ (abs, SF.polygamma, SF.trigamma, acsc, acsch, asech, NaNMath.log2, NaNMath.log10)
        continue
    end

    fun = eval(:(ξ ->($f)(ξ)))

    fd = ForwardDiff.derivative(fun, x)
    sym = Symbolics.Differential(x)(fun(x)) |> expand_derivatives

    @test isequal(fd, sym)
end

# These are evaluated numerically. For some reason isequal evaluates to false for the symbolic expressions.
for f ∈ (acsc, asech, NaNMath.log2, NaNMath.log10)
    fun = eval(:(ξ ->($f)(ξ)))

    fd = ForwardDiff.derivative(fun, 1.0)
    sym = Symbolics.Differential(x)(fun(x)) |> expand_derivatives

    @test fd ≈ substitute(sym, Dict(x => 1.0))
end

for f ∈ SymbolicUtils.basic_diadic
    if f ∈ (//,)
        continue
    end

    fun = eval(:(ξ ->($f)(ξ, 2.0)))

    fd = ForwardDiff.derivative(fun, x)
    sym = Symbolics.Differential(x)(fun(x)) |> expand_derivatives

    @test isequal(fd, sym)
end

for f ∈ SymbolicUtils.diadic
    if f ∈ (max, min, NaNMath.atanh, mod, rem, copysign, besselj, bessely, besseli, besselk)
        continue
    end

    fun = eval(:(ξ ->($f)(ξ, 2.0)))

    fd = ForwardDiff.derivative(fun, x)
    sym = Symbolics.Differential(x)(fun(x)) |> expand_derivatives

    @test isequal(fd, sym)
end

for f ∈ (NaNMath.atanh,)
    fun = eval(:(ξ ->($f)(ξ)))

    fd = ForwardDiff.derivative(fun, x)
    sym = Symbolics.Differential(x)(fun(x)) |> expand_derivatives

    @test isequal(fd, sym)
end

for f ∈ (besselj, bessely, besseli, besselk)
    fun = eval(:(ξ ->($f)(ξ, 2)))

    fd = ForwardDiff.derivative(fun, x)
    sym = Symbolics.Differential(x)(fun(x)) |> expand_derivatives

    @test isequal(fd, sym)
end

# Additionally test these definitions from ForwardDiff #
#------------------------------------------------------#

# https://github.com/JuliaDiff/ForwardDiff.jl/blob/d3002093beb88ff0b98ed178377961dfd55c1247/src/dual.jl#L599
# and
# https://github.com/JuliaDiff/ForwardDiff.jl/blob/d3002093beb88ff0b98ed178377961dfd55c1247/src/dual.jl#L683
for f ∈ (hypot, muladd)
    fun = eval(:(ξ ->($f)(ξ, 2.0, 3.0)))

    fd = ForwardDiff.derivative(fun, 5.0)
    sym = Symbolics.Differential(x)(fun(x)) |> expand_derivatives

    @test fd ≈ substitute(sym, Dict(x => 5.0))
end

# fma is not defined for Symbolics.Num

# https://github.com/JuliaSymbolics/Symbolics.jl/issues/1246
@testset "isequal type ambiguity" begin
    @variables z
    y(x) = isequal(z, x) ? 0 : x
    @test ForwardDiff.derivative(y, 0) == 1 # expect ∂(x)/∂x
end

@testset "NaNMath.pow (issue #1399)" begin
    @variables x
    @test_throws DomainError substitute(ForwardDiff.derivative(z -> x^z, 0.5), x => -1.0)
    @test isnan(Symbolics.value(substitute(ForwardDiff.derivative(z -> NaNMath.pow(x, z), 0.5), x => -1.0)))
end
