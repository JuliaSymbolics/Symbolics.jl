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

struct RLArrField
    v::Vector{Float64}
end
@symstruct RLArrField

# Nine fields, one more than `RECORD_LITERAL_MAX_FIELDS`, so no methods are generated.
struct RLWide
    a1::Int; a2::Int; a3::Int; a4::Int; a5::Int
    a6::Int; a7::Int; a8::Int; a9::Int
end
@symstruct RLWide

# Opts out of the generated constructors and handles the symbolic case itself, which is
# what a struct's owner can do to avoid `2^nfields - 1` methods.
struct RLHand
    a::Float64
    b::Float64
    function RLHand(args...)
        any(x -> unwrap(x) isa Symbolics.SymbolicT, args) &&
            return record_literal(RLHand, args)
        return new(args...)
    end
end
@symstruct RLHand begin
    literal_constructors() = false
end

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
    @test is_record_literal(record_literal(RLWide, args))
    @test_throws ArgumentError record_literal(RLPair, (1.0,))
end

@testset "applying the operation does not rely on the generated constructors" begin
    # Concrete fields give the value itself, for a struct of either width.
    @test Symbolics.RecordLiteral{RLPair}()(1.0, 2.0) == RLPair(1.0, 2.0)
    @test Symbolics.RecordLiteral{RLWide}()(ntuple(i -> i, 9)...) == RLWide(ntuple(i -> i, 9)...)

    # A symbolic field gives a literal. For `RLPair` the registered constructor would do
    # this anyway; `RLWide` has none, so the operation has to build it itself.
    @test is_record_literal(unwrap(Symbolics.RecordLiteral{RLPair}()(unwrap(q1), 2.0)))
    wide = Symbolics.RecordLiteral{RLWide}()(unwrap(q1), ntuple(i -> i + 1, 8)...)
    @test is_record_literal(unwrap(wide))
    @test SU.symtype(unwrap(wide)) === RLWide
end

@testset "a single array-valued field is not read as a field list" begin
    # The field list is a `Tuple`; anything else is one field value. Without that, an
    # array-valued field of a one-field struct is mistaken for two fields.
    lit = unwrap(record_literal(RLArrField, [1.0, 2.0]))
    @test is_record_literal(lit)
    @test length(SU.arguments(lit)) == 1
    @test SU.symtype(lit) === RLArrField

    # A tuple is still read as the field list.
    @test length(SU.arguments(unwrap(record_literal(RLPair, (1.0, 2.0))))) == 2
end

@testset "`literal_constructors() = false` skips the generated methods" begin
    # Registration is unaffected; only the constructor methods are skipped.
    @test Symbolics.has_symwrapper(RLHand)
    @test length(methods(RLHand)) == 1
    @test length(methods(RLPair)) > 1

    # The hand-written constructor covers both paths itself.
    @test RLHand(1.0, 2.0) isa RLHand
    lit = RLHand(q1, 2.0)
    @test is_record_literal(unwrap(lit))
    @test SU.symtype(unwrap(lit)) === RLHand
    @test isequal(unwrap(SymStruct{RLHand}(unwrap(lit)).a), unwrap(q1))
end

@testset "field symtypes are checked against the field types" begin
    # A symbolic is usually declared more loosely than the field it fills, so this has to
    # be accepted: `@variables` gives symtype `Real` for a `Float64` field.
    @test is_record_literal(RLPair(q1, q2))
    @test is_record_literal(RLPair(1.0, q2))

    # An unrelated symtype is a mistake in either direction. Caught here rather than on
    # field access, where it used to surface as an assertion inside `getproperty`.
    @variables rec::RLPair str::String
    @test_throws ArgumentError RLPair(rec, q1)
    @test_throws ArgumentError RLPair(str, q1)
    @test_throws ArgumentError record_literal(RLPair, (rec, q1))

    # The message names the field, its type, and the symtype received.
    err = try
        RLPair(rec, q1)
    catch e
        sprint(showerror, e)
    end
    @test occursin("field `x`", err)
    @test occursin("RLPair", err)
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
