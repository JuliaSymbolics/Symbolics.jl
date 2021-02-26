module Symbolics

using DocStringExtensions

using LinearAlgebra

using Reexport

@reexport using SymbolicUtils

import SymbolicUtils: Term, Add, Mul, Pow, Sym, symtype,
                      FnType, @rule, Rewriters, substitute, similarterm,
                      promote_symtype, istree, operation, arguments

import SymbolicUtils.Rewriters: Chain, Prewalk, Postwalk, Fixpoint

# re-export

export simplify, substitute

using SciMLBase
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
