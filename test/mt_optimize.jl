using Symbolics
using SymbolicUtils
using Test 

@syms a b

# cmp(x,y) = x == y

# function cmp(x::SymbolicUtils.Term, y::SymbolicUtils.Term)
#     lx, ly = length(x.arguments), length(y.arguments) 
#     return x.f == y.f && lx == ly &&
#         all([cmp(x.arguments[i], y.arguments[i]) for i in 1:lx])
# end

ex = 2a + 2b - (a*(a + b))
res = Symbolics.optimize(ex)
