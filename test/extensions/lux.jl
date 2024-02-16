using Symbolics, Lux, Random, Test
using ComponentArrays

@testset "Dense" begin
    model = Dense(5, 6)
    rng = Random.default_rng()
    x = randn(rng, Float32, 5)
    ps, st = LuxCore.setup(rng, model)

    ca = ComponentArray(ps)
    Symbolics.@variables sym_ca[1:length(ca)] = ca
    Symbolics.@variables sym_ps::typeof(ps) = ps
    Symbolics.@variables sym_st::typeof(st) = st
    Symbolics.@variables sym_x[1:5] = Float32[1,2,3,4,5]

    out_ref = LuxCore.partial_apply(model, x, ps, st)
    @test out_ref isa Vector{Float32}

    out = LuxCore.partial_apply(model, sym_x, ps, st)
    # if the symbolic function regstration wouldn't work, we'd get the
    # vector propagated through the neural network and get a Vector{Num}
    @test out isa Symbolics.Arr
    @test length(out) == 6
    # test that we  can recover the same value as when using concrete numbers
    out_sub = Symbolics.value.(Symbolics.substitute.(Symbolics.scalarize(out), (Dict(sym_x => x),)))
    @test out_sub == out_ref

    out = LuxCore.partial_apply(model, sym_x, sym_ps, st)
    @test out isa Symbolics.Arr
    @test length(out) == 6
    out_sub = Symbolics.value.(Symbolics.substitute.(Symbolics.scalarize(out), (Dict(sym_x => x, sym_ps=>ps),)))
    @test out_sub == out_ref

    out = LuxCore.partial_apply(model, sym_x, ps, sym_st)
    @test out isa Symbolics.Arr
    @test length(out) == 6
    out_sub = Symbolics.value.(Symbolics.substitute.(Symbolics.scalarize(out), (Dict(sym_x => x, sym_st => st),)))
    @test out_sub == out_ref

    out = LuxCore.partial_apply(model, sym_x, sym_ps, sym_st)
    @test out isa Symbolics.Arr
    @test length(out) == 6
    out_sub = Symbolics.value.(Symbolics.substitute.(Symbolics.scalarize(out), (Dict(sym_x => x, sym_ps => ps, sym_st => st),)))
    @test out_sub == out_ref

    out = LuxCore.partial_apply(model, sym_x, ca, st)
    @test out isa Symbolics.Arr
    @test length(out) == 6
    out_sub = Symbolics.value.(Symbolics.substitute.(Symbolics.scalarize(out), (Dict(sym_x => x,),)))
    @test out_sub == out_ref

    out = LuxCore.partial_apply(model, sym_x, sym_ca, st)
    @test out isa Symbolics.Arr
    @test length(out) == 6
    out_sub = Symbolics.value.(Symbolics.substitute.(Symbolics.scalarize(out), (Dict(sym_x => x, sym_ca => ca),)))
    @test out_sub == out_ref

    out = LuxCore.partial_apply(model, sym_x, sym_ca, sym_st)
    @test out isa Symbolics.Arr
    @test length(out) == 6
    out_sub = Symbolics.value.(Symbolics.substitute.(Symbolics.scalarize(out), (Dict(sym_x => x, sym_ca => ca, sym_st => st),)))
    @test out_sub == out_ref
end

@testset "Chain" begin
    model = Chain(Dense(5, 6), Dense(6, 2), Dense(2, 3))
    rng = Random.default_rng()
    x = randn(rng, Float32, 5)
    ps, st = LuxCore.setup(rng, model)

    ca = ComponentArray(ps)
    Symbolics.@variables sym_ca[1:length(ca)] = ca
    Symbolics.@variables sym_ps::typeof(ps) = ps
    Symbolics.@variables sym_st::typeof(st) = st
    Symbolics.@variables sym_x[1:5] = Float32[1, 2, 3, 4, 5]

    out_ref = LuxCore.partial_apply(model, x, ps, st)
    @test out_ref isa Vector{Float32}

    out = LuxCore.partial_apply(model, sym_x, ps, st)
    # if the symbolic function regstration wouldn't work, we'd get the
    # vector propagated through the neural network and get a Vector{Num}
    @test out isa Symbolics.Arr
    @test length(out) == 3
    # test that we  can recover the same value as when using concrete numbers
    out_sub = Symbolics.value.(Symbolics.substitute.(Symbolics.scalarize(out), (Dict(sym_x => x),)))
    @test out_sub == out_ref

    out = LuxCore.partial_apply(model, sym_x, sym_ps, st)
    @test out isa Symbolics.Arr
    @test length(out) == 3
    out_sub = Symbolics.value.(Symbolics.substitute.(Symbolics.scalarize(out), (Dict(sym_x => x, sym_ps => ps),)))
    @test out_sub == out_ref

    out = LuxCore.partial_apply(model, sym_x, ps, sym_st)
    @test out isa Symbolics.Arr
    @test length(out) == 3
    out_sub = Symbolics.value.(Symbolics.substitute.(Symbolics.scalarize(out), (Dict(sym_x => x, sym_st => st),)))
    @test out_sub == out_ref

    out = LuxCore.partial_apply(model, sym_x, sym_ps, sym_st)
    @test out isa Symbolics.Arr
    @test length(out) == 3
    out_sub = Symbolics.value.(Symbolics.substitute.(Symbolics.scalarize(out), (Dict(sym_x => x, sym_ps => ps, sym_st => st),)))
    @test out_sub == out_ref

    out = LuxCore.partial_apply(model, sym_x, ca, st)
    @test out isa Symbolics.Arr
    @test length(out) == 3
    out_sub = Symbolics.value.(Symbolics.substitute.(Symbolics.scalarize(out), (Dict(sym_x => x,),)))
    @test out_sub == out_ref

    out = LuxCore.partial_apply(model, sym_x, sym_ca, st)
    @test out isa Symbolics.Arr
    @test length(out) == 3
    out_sub = Symbolics.value.(Symbolics.substitute.(Symbolics.scalarize(out), (Dict(sym_x => x, sym_ca => ca),)))
    @test out_sub == out_ref

    out = LuxCore.partial_apply(model, sym_x, sym_ca, sym_st)
    @test out isa Symbolics.Arr
    @test length(out) == 3
    out_sub = Symbolics.value.(Symbolics.substitute.(Symbolics.scalarize(out), (Dict(sym_x => x, sym_ca => ca, sym_st => st),)))
    @test out_sub == out_ref
end
