module Symbolics

using DocStringExtensions

using LinearAlgebra

import SymbolicUtils

import SymbolicUtils: Term, Add, Mul, Pow, Sym, symtype, to_symbolic,
                      FnType, @rule, Rewriters, substitute, similarterm,
                      promote_symtype, istree, operation, arguments


# re-export

export simplify, substitute

include("utils.jl")

export Num
include("num.jl")

using TreeViews
export @variables
include("variable.jl")

include("linearity.jl")

using DiffRules, SpecialFunctions, NaNMath
export Differential, expand_derivatives, gradient,
       jacobian, jacobian_sparsity, sparsejacobian,
       hessian, sparsehessian, hessian_sparsity

include("diff.jl")

end # module
