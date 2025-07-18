using SnoopCompile: @snoop_invalidations
using Symbolics

struct FakeType end

@testset "`option_to_metadata_type`" begin
    invs = @snoop_invalidations Symbolics.option_to_metadata_type(::Val{:____!_internal_999}) = FakeType
    @test isempty(invs)
end
