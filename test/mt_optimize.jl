using Symbolics
using SymbolicUtils
using Test 

@syms a b
tr(f,x,y) = SymbolicUtils.Term{Real}(f, [x,y])
tn(f,x,y) = SymbolicUtils.Term{Number}(f, [x,y])


ex = 2a + 2b - (a*(a + b))
prex = Symbolics.preprocess(ex)
res = Symbolics.optimize(ex)

@test isequal(res, tn(*, tn(+,b,a), tn(-, 2, a)))


res = Symbolics.optimize(sin(a^2)/cos(a^2))
@test isequal(res, tan(a^2))