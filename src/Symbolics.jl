module Symbolics

using DocStringExtensions

using LinearAlgebra

import SymbolicUtils

import SymbolicUtils: Term, Add, Mul, Pow, Sym, symtype, to_symbolic,
                      FnType, @rule, Rewriters, substitute, similarterm,
                      promote_symtype, istree, operation, arguments

import SymbolicUtils.Rewriters: Chain, Prewalk, Postwalk, Fixpoint

# re-export

export simplify, substitute

export Num
include("num.jl")

include("utils.jl")

using MacroTools
include("register.jl")

using TreeViews
export @variables
include("variable.jl")

include("linearity.jl")

using DiffRules, SpecialFunctions, NaNMath

using SparseArrays

export Differential, expand_derivatives, gradient,
       jacobian, jacobian_sparsity, sparsejacobian,
       hessian, sparsehessian, hessian_sparsity

include("diff.jl")

include("linear_algebra.jl")

end # module
