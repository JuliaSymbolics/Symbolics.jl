using Symbolics
using ADTypes: ADTypes
using Test

detector = SymbolicsSparsityDetector()
@test detector isa ADTypes.AbstractSparsityDetector

f(x) = reshape(vcat(x[1], diff(vec(x))), size(x))
f!(y, x) = copyto!(vec(y), vec(f(x)))

for x in (rand(4), rand(4, 5))
    @test sum(ADTypes.jacobian_sparsity(f, x, detector)) == 2length(x) - 1
    @test sum(ADTypes.jacobian_sparsity(f!, similar(x), x, detector)) == 2length(x) - 1
end

g(x) = sum(abs2, diff(vec(x)))

for x in (rand(4), rand(4, 5))
    @test sum(ADTypes.hessian_sparsity(g, x, detector)) == 3length(x) - 2
end
