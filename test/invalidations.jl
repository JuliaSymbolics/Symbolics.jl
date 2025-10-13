using SnoopCompile: @snoop_invalidations, invalidation_trees
using Symbolics

struct FakeType end

@testset "`option_to_metadata_type`" begin
    invs = @snoop_invalidations Symbolics.option_to_metadata_type(::Val{:____!_internal_999}) = FakeType
    @test isempty(invalidation_trees(invs))
end
