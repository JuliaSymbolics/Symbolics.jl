using SciMLTesting
import StaticArrays
using Test

symbolic_utils = Symbolics.SymbolicUtils
basic_symbolic = symbolic_utils.BasicSymbolic{symbolic_utils.SymReal}
basic_symbolic_wrapper = getfield(parentmodule(basic_symbolic), nameof(basic_symbolic))

# Upstream compatibility names that are not declared public by their owners.
# Issue #2014 names (`map`, `parse`, `require_one_based_indexing`, `reverse`) must
# NOT be added here — fix call sites to use public APIs or local helpers instead.
const QUALIFIED_ACCESS_IGNORE = (
    :BlasInt, :Cartesian, :Experimental, :ParseError, :ReshapedArray,
    :Unknown, :acos, :acosh, :alignment, :asin, :atanh, :checknonsingular,
    :cos, :diffrule, :diffrules, :eval, :getdoc, :log, :log10, :log1p,
    :log2, :max, :min, :nocolor, :power_by_squaring, :register_error_hint,
    :AbstractCompressedVector, :AbstractSparseMatrixCSC, :AbstractTriangular,
    :Slice, :StaticArray, :TwicePrecision, :TypedEndpointsInterval,
    :sin, :sqrt, :striplines, :tan,
)

# Names from https://github.com/JuliaSymbolics/Symbolics.jl/issues/2014
const ISSUE_2014_NAMES = (:map, :parse, :require_one_based_indexing, :reverse)

@testset "issue #2014: stdlib accesses use public APIs" begin
    for name in ISSUE_2014_NAMES
        @test name ∉ QUALIFIED_ACCESS_IGNORE
    end
    # On Julia 1.11 these names are not public; call sites must not qualify them.
    # (On 1.12+ they are public, so this filter is vacuously empty.)
    # ExplicitImports is a transitive test dep via SciMLTesting; load it from
    # `Base.require` so the test env does not need a direct entry.
    if VERSION >= v"1.11"
        ei_id = Base.PkgId(
            Base.UUID("7d51a73a-1435-4ff3-83d9-f097790105c7"),
            "ExplicitImports",
        )
        ExplicitImports = Base.require(ei_id)
        offending = Symbol[]
        for (_, rows) in ExplicitImports.improper_qualified_accesses(Symbolics)
            for row in rows
                if !row.public_access && row.name in ISSUE_2014_NAMES
                    push!(offending, row.name)
                end
            end
        end
        @test isempty(offending)
    end
end

run_qa(
    Symbolics;
    aqua_kwargs = (;
        # Shared SymbolicUtils interfaces are jointly owned. Base has no binary `~`, so
        # Symbolics' equation operator cannot collide with anything it does not own.
        piracies = (;
            treat_as_own = (
                Base.:~,
                basic_symbolic_wrapper,
                symbolic_utils.arguments,
                symbolic_utils.Code.cse_inside_expr,
                symbolic_utils.promote_shape,
                symbolic_utils.promote_symtype,
            ),
        ),
    ),
    ei_kwargs = (;
        # These are upstream names used for compatibility with Base, LinearAlgebra,
        # DiffRules, MacroTools, and NaNMath; those owners do not declare them public.
        all_qualified_accesses_are_public = (;
            ignore = QUALIFIED_ACCESS_IGNORE,
        ),
    ),
    reexports_allow = (
        Symbol("@acrule"), Symbol("@arrayop"), Symbol("@makearray"),
        Symbol("@rule"), Symbol("@syms"), :BS, :IRStructure, :Rewriters,
        :RuleSet, :SafeReal, :SymReal, :SymbolicUtils, :TreeReal, :Unknown,
        :arguments, :expand, :flatten_fractions, :get_reachability, :getmetadata,
        :hasmetadata, :ifelse_branching, :ifelse_eager, :iscall, :istree,
        :operation, :populate_ir!, :print_ir, :quick_cancel, :scalarize,
        :setmetadata, :shape, :simplify, :simplify_fractions, :sorted_arguments,
        :substitute, :term, :unwrap, :unwrap_const, :vartype,
    ),
)
