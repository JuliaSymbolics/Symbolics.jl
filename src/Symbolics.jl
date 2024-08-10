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
import TermInterface: maketerm, iscall, operation, arguments, symtype, metadata

import SymbolicUtils: Term, Add, Mul, Pow, Sym, Div, BasicSymbolic,
FnType, @rule, Rewriters, substitute,
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
include("utils.jl")

using ConstructionBase
include("arrays.jl")

export @register_symbolic, @register_array_symbolic
include("register.jl")

export @variables, Variable
include("variable.jl")

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

include("solver.jl")
export solve_single_eq
export solve_system_eq
export lambertw

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
include("./solver/coeffs.jl")
include("./solver/nemo_stuff.jl")
include("./solver/solve_helpers.jl")
include("./solver/postprocess.jl")
include("./solver/univar.jl")
include("./solver/ia_helpers.jl")
include("./solver/polynomialization.jl")
include("./solver/attract.jl")
include("./solver/ia_main.jl")
include("./solver/main.jl")

function symbolics_to_sympy end
export symbolics_to_sympy

include("../ext/SymbolicsForwardDiffExt.jl")
using ..SymbolicsForwardDiffExt

@static if !isdefined(Base, :get_extension)
    using Requires
end

@static if !isdefined(Base,:get_extension)
    function __init__()
        @require Groebner="0b43b601-686d-58a3-8a1c-6623616c7cd4" begin
            include("../ext/SymbolicsGroebnerExt.jl")
        end
        @require SymPy="24249f21-da20-56a4-8eb1-6a02cf4ae2e6" begin
            include("../ext/SymbolicsSymPyExt.jl")
        end
        @require Nemo="2edaba10-b0f1-5616-af89-8c11ac63239a" begin
            include("../ext/SymbolicsNemoExt.jl")
        end
    end
end

end # module
