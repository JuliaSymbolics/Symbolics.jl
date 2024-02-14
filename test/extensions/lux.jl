using Symbolics, Lux, Random, Test

@testset "Symbolics extension" begin
    model = Dense(5, 6)
    rng = Random.default_rng()
    x = randn(rng, Float32, 5)
    ps, st = LuxCore.setup(rng, model)

    Symbolics.@variables sym_ps[1:5] = Float32[1, 2, 3, 4, 5]

    out = LuxCore.partial_apply(model, sym_ps, ps, st)
    @test out isa Symbolics.Arr
end
