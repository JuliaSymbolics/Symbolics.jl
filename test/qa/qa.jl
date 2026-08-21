using SciMLTesting
using Symbolics

run_qa(
    Symbolics;
    ei_kwargs = (;
        # These are upstream names used for compatibility with Base, LinearAlgebra,
        # DiffRules, MacroTools, and NaNMath; those owners do not declare them public.
        all_qualified_accesses_are_public = (; ignore = (
            :BlasInt, :Cartesian, :Experimental, :ParseError, :ReshapedArray,
            :Unknown, :acos, :acosh, :alignment, :asin, :atanh, :checknonsingular,
            :cos, :diffrule, :diffrules, :eval, :getdoc, :log, :log10, :log1p,
            :log2, :max, :min, :nocolor, :power_by_squaring, :register_error_hint,
            :sin, :sqrt, :striplines, :tan,
        )),
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
