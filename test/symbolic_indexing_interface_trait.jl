using Symbolics
using SymbolicUtils
using SymbolicIndexingInterface

@test all(symbolic_type.([SymbolicUtils.BasicSymbolic, Symbolics.Num]) .==
          (ScalarSymbolic(),))
@test symbolic_type(Symbolics.Arr) == ArraySymbolic()
@variables x
@test symbolic_type(x) == ScalarSymbolic()
@variables y[1:3]
@test symbolic_type(y) == ArraySymbolic()
@test all(symbolic_type.(collect(y)) .== (ScalarSymbolic(),))
