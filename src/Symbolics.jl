module Symbolics

import SymbolicUtils

import SymbolicUtils: Term, Add, Mul, Pow, Sym, symtype, to_symbolic,
                      FnType, @rule, Rewriters, substitute, similarterm,
                      promote_symtype

export Num

using LinearAlgebra

include("utils.jl")

include("num.jl")

export @variables
include("variable.jl")

end # module
