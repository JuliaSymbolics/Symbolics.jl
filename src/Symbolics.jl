"""
$(DocStringExtensions.README)
"""
module Symbolics

using PrecompileTools

import PrecompileTools: @recompile_invalidations

@recompile_invalidations begin
    import CommonWorldInvalidations
end

using DocStringExtensions, Markdown

using LinearAlgebra

using Primes

using Reexport

using DomainSets

using Setfield

import DomainSets: Domain

using TermInterface
import TermInterface: maketerm, iscall, operation, arguments, metadata

import SymbolicUtils: Term, Add, Mul, Pow, Sym, Div, BasicSymbolic,
FnType, @rule, Rewriters, substitute, symtype,
promote_symtype, isadd, ismul, ispow, isterm, issym, isdiv

using SymbolicUtils.Code

import SymbolicUtils.Rewriters: Chain, Prewalk, Postwalk, Fixpoint

import SymbolicUtils.Code: toexpr

import ArrayInterface
using RuntimeGeneratedFunctions
using SciMLBase, IfElse
import MacroTools

using SymbolicIndexingInterface

import SymbolicLimits

using ADTypes: ADTypes

@reexport using SymbolicUtils
RuntimeGeneratedFunctions.init(@__MODULE__)

# re-export

export simplify, substitute

export Num
import MacroTools: splitdef, combinedef, postwalk, striplines
include("wrapper-types.jl")

include("num.jl")

include("rewrite-helpers.jl")
include("complex.jl")

"""
    substitute(expr, s)

Performs the substitution on `expr` according to rule(s) `s`.
# Examples
```julia
julia> @variables t x y z(t)
4-element Vector{Num}:
    t
    x
    y
 z(t)
julia> ex = x + y + sin(z)
(x + y) + sin(z(t))
julia> substitute(ex, Dict([x => z, sin(z) => z^2]))
(z(t) + y) + (z(t) ^ 2)
```
"""
substitute

export Equation, ConstrainedEquation
include("equations.jl")

export Inequality, ≲, ≳
include("inequality.jl")

import Bijections, DynamicPolynomials
export tosymbol
include("utils.jl")

using ConstructionBase
include("arrays.jl")

export @register_symbolic, @register_array_symbolic
include("register.jl")

export @variables, Variable
include("variable.jl")

function slog end; function ssqrt end; function scbrt end
include("linearity.jl")

using DiffRules, SpecialFunctions, NaNMath

using SparseArrays

export Differential, expand_derivatives, is_derivative

include("diff.jl")

export SymbolicsSparsityDetector

include("adtypes.jl")

export Difference, DiscreteUpdate

include("difference.jl")

export infimum, supremum
include("domains.jl")

export Integral
include("integral.jl")

include("array-lib.jl")

using LogExpFunctions
include("logexpfunctions-lib.jl")

include("linear_algebra.jl")
export symbolic_linear_solve, solve_for

include("groebner_basis.jl")
export groebner_basis, is_groebner_basis

import Libdl
include("build_function.jl")
export build_function

import Distributions
include("extra_functions.jl")

using Latexify
using LaTeXStrings
include("latexify_recipes.jl")

using RecipesBase
include("plot_recipes.jl")

include("semipoly.jl")


include("parsing.jl")
export parse_expr_to_symbolic

include("error_hints.jl")
include("struct.jl")
include("operators.jl")

include("limits.jl")
export limit

# Hacks to make wrappers "nicer"
const NumberTypes = Union{AbstractFloat,Integer,Complex{<:AbstractFloat},Complex{<:Integer}}
(::Type{T})(x::SymbolicUtils.Symbolic) where {T<:NumberTypes} = throw(ArgumentError("Cannot convert Sym to $T since Sym is symbolic and $T is concrete. Use `substitute` to replace the symbolic unwraps."))
for T in [Num, Complex{Num}]
    @eval begin
        #(::Type{S})(x::$T) where {S<:Union{NumberTypes,AbstractArray}} = S(Symbolics.unwrap(x))::S

        SymbolicUtils.simplify(n::$T; kw...) = wrap(SymbolicUtils.simplify(unwrap(n); kw...))
        SymbolicUtils.simplify_fractions(n::$T; kw...) = wrap(SymbolicUtils.simplify_fractions(unwrap(n); kw...))
        SymbolicUtils.expand(n::$T) = wrap(SymbolicUtils.expand(unwrap(n)))
        substitute(expr::$T, s::Pair; kw...) = wrap(substituter(s)(unwrap(expr); kw...)) # backward compat
        substitute(expr::$T, s::Vector; kw...) = wrap(substituter(s)(unwrap(expr); kw...))
        substitute(expr::$T, s::Dict; kw...) = wrap(substituter(s)(unwrap(expr); kw...))

        SymbolicUtils.Code.toexpr(x::$T) = SymbolicUtils.Code.toexpr(unwrap(x))

        SymbolicUtils.setmetadata(x::$T, t, v) = wrap(SymbolicUtils.setmetadata(unwrap(x), t, v))
        SymbolicUtils.getmetadata(x::$T, t) = SymbolicUtils.getmetadata(unwrap(x), t)
        SymbolicUtils.hasmetadata(x::$T, t) = SymbolicUtils.hasmetadata(unwrap(x), t)

        Broadcast.broadcastable(x::$T) = x
    end
    for S in [:(Symbolic{<:FnType}), :CallWithMetadata]
        @eval (f::$S)(x::$T, y...) = wrap(f(unwrap(x), unwrap.(y)...))
    end
end

for sType in [Pair, Vector, Dict]
    @eval substitute(expr::Arr, s::$sType; kw...) = wrap(substituter(s)(unwrap(expr); kw...))
end

# Symbolic solver
include("solver/preprocess.jl")
include("solver/nemo_stuff.jl")
include("solver/solve_helpers.jl")
include("solver/postprocess.jl")
include("solver/univar.jl")
include("solver/ia_helpers.jl")
include("solver/polynomialization.jl")
include("solver/attract.jl")
include("solver/ia_main.jl")
include("solver/main.jl")
include("solver/new_feature.jl")
export symbolic_solve

function symbolics_to_sympy end
export symbolics_to_sympy

function __init__()
    Base.Experimental.register_error_hint(TypeError) do io, exc
        if exc.expected == Bool && exc.got isa Num
            println(io,
                "\nA symbolic expression appeared in a Boolean context. This error arises in situations where Julia expects a Bool, like ")
            printstyled(io, "if boolean_condition", color = :blue)
            printstyled(
                io, "\t\t use ifelse(boolean_condition, then branch, else branch)\n",
                color = :green)
            printstyled(io, "x && y", color = :blue)
            printstyled(io, "\t\t\t\t use x & y\n", color = :green)
            printstyled(io, "boolean_condition ? a : b", color = :blue)
            printstyled(io, "\t use ifelse(boolean_condition, a, b)\n", color = :green)
            print(io,
                "but a symbolic expression appeared instead of a Bool. For help regarding control flow with symbolic variables, see https://docs.sciml.ai/ModelingToolkit/dev/basics/FAQ/#How-do-I-handle-if-statements-in-my-symbolic-forms?")
        end
    end
end

end # module
