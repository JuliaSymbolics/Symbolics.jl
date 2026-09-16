using SciMLTesting
import StaticArrays

symbolic_utils = Symbolics.SymbolicUtils
basic_symbolic = symbolic_utils.BasicSymbolic{symbolic_utils.SymReal}
basic_symbolic_wrapper = getfield(parentmodule(basic_symbolic), nameof(basic_symbolic))

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
        # `DefaultSubstituter` is the concrete substitution protocol needed to let a
        # `Num` widen after complex replacement; `promote_op(matprod, ...)` is Julia's
        # element-type inference hook used to keep SymbolicNumber matrix products concrete.
        all_qualified_accesses_are_public = (;
            ignore = (
                :BlasInt, :Cartesian, :DefaultSubstituter, :Experimental, :ParseError,
                :ReshapedArray, :Unknown, :acos, :acosh, :alignment, :asin, :atanh,
                :checknonsingular, :cos, :diffrule, :diffrules, :eval, :getdoc, :log,
                :log10, :log1p, :log2, :matprod, :max, :min, :nocolor,
                :power_by_squaring, :promote_op, :register_error_hint,
                :AbstractCompressedVector, :AbstractSparseMatrixCSC, :AbstractTriangular,
                :Slice, :StaticArray, :TwicePrecision, :TypedEndpointsInterval,
                :sin, :sqrt, :striplines, :tan,
            ),
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
