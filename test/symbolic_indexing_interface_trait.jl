using Symbolics
using SymbolicUtils
using SymbolicIndexingInterface

@test all(symbolic_type.([SymbolicUtils.BasicSymbolic, Symbolics.Num, Symbolics.CallWithMetadata]) .==
          (ScalarSymbolic(),))
@test all(symbolic_type.([Symbolics.ArrayOp, Symbolics.Arr]) .==
          (ArraySymbolic(),))
@variables x
@test symbolic_type(x) == ScalarSymbolic()
@variables y[1:3]
@test symbolic_type(y) == ArraySymbolic()
@test all(symbolic_type.(collect(y)) .== (ScalarSymbolic(),))
@test symbolic_type(Symbolics.unwrap(y .* y)) == ArraySymbolic()
@variables z(..)
@test symbolic_type(z) == ScalarSymbolic()

@variables x y z
subs = Dict(x => 0.1, y => 2z)
subs2 = merge(subs, Dict(z => 2x+3))

@test symbolic_evaluate(x, subs) == 0.1
@test isequal(symbolic_evaluate(y, subs), 2z)
@test symbolic_evaluate(y, subs2) == 6.4
