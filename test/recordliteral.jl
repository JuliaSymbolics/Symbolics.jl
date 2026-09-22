using Symbolics
using SymbolicUtils, SymbolicUtils.Code
import SymbolicUtils as SU
using Symbolics: unwrap, is_record_literal, record_literal
using Test

struct RLPair
    x::Float64
    y::Float64
end
@symstruct RLPair

struct RLNested
    inner::RLPair
    k::Float64
end
@symstruct RLNested

struct RLValidated
    a::Float64
    RLValidated(a) = (a < 0 && error("a must be >= 0"); new(a))
end
@symstruct RLValidated

struct RLParam{V}
    x::Int
    z::V
end
@symstruct RLParam{V}

# Nine fields, one more than `RECORD_LITERAL_MAX_FIELDS`, so no methods are generated.
struct RLWide
    a1::Int; a2::Int; a3::Int; a4::Int; a5::Int
    a6::Int; a7::Int; a8::Int; a9::Int
end
@symstruct RLWide

@variables q1 q2 q3

@testset "construction" begin
    # Fully concrete construction is untouched.
    @test RLPair(1.0, 2.0) isa RLPair

    lit = unwrap(RLPair(q1, q2))
    @test is_record_literal(lit)
    @test SU.symtype(lit) === RLPair
    @test SU.shape(lit) == SU.ShapeVecT()
    @test length(SU.arguments(lit)) == 2

    # A single symbolic argument is enough.
    @test is_record_literal(unwrap(RLPair(1.0, q2)))
    @test is_record_literal(unwrap(RLPair(q1, 2.0)))

    @test SU.promote_symtype(Symbolics.RecordLiteral{RLPair}(), Float64, Float64) === RLPair
end

@testset "the struct's own constructor is preserved" begin
    @test RLValidated(1.0) isa RLValidated
    # Validation still runs for concrete construction ...
    @test_throws ErrorException RLValidated(-1.0)
    # ... and is skipped for a literal, which never calls the constructor.
    @test is_record_literal(unwrap(RLValidated(q1)))
end

@testset "parametric structs" begin
    @test RLParam{Float64}(1, 2.0) isa RLParam{Float64}
    lit = unwrap(RLParam{Float64}(q1, 2.0))
    @test is_record_literal(lit)
    @test SU.symtype(lit) === RLParam{Float64}
end

@testset "wide structs fall back to `record_literal`" begin
    args = ntuple(i -> i, 9)
    # No constructor methods are generated, so this stays concrete.
    @test RLWide(args...) isa RLWide
    @test is_record_literal(unwrap(record_literal(RLWide, args)))
    @test_throws ArgumentError record_literal(RLPair, (1.0,))
end

@testset "field access folds through a literal" begin
    lit = RLPair(q1, q2)
    @test isequal(unwrap(SymStruct{RLPair}(unwrap(lit)).x), unwrap(q1))
    @test isequal(unwrap(SymStruct{RLPair}(unwrap(lit)).y), unwrap(q2))

    nlit = RLNested(lit, q3)
    wrapped = SymStruct{RLNested}(unwrap(nlit))
    @test isequal(unwrap(wrapped.k), unwrap(q3))
    # Folding composes: the inner literal folds in turn.
    @test isequal(unwrap(wrapped.inner.x), unwrap(q1))
end

@testset "show" begin
    @test occursin("RLPair(", sprint(show, unwrap(RLPair(q1, q2))))
end

@testset "code generation" begin
    # A literal lowers to a call to the struct's constructor.
    f = eval(build_function(RLPair(q1, q2), q1, q2; expression = Val{true}))
    @test Base.invokelatest(f, 1.0, 2.0) == RLPair(1.0, 2.0)

    # A field projected out of a literal generates that field's expression.
    g = eval(build_function(
        SymStruct{RLPair}(unwrap(RLPair(q1 + 1, q2))).x, q1, q2; expression = Val{true}))
    @test Base.invokelatest(g, 1.0, 2.0) == 2.0
end
