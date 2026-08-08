using SymbolicUtils.Rewriters
using LinearAlgebra

isarray(x) = symtype(x) <: AbstractArray
ismatrix(x) = ndims(x) == 2
isvector(x) = ndims(x) == 1

array_rules = [@rule adjoint(adjoint(~X)) => ~X
               @rule ((~A*~B)*~C) => (~A*(~B*~C)) where size(~A,1)*size(~B,2) > size(~B,1)*size(~C,2)
               @rule broadcast(~f, broadcast(~g, ~A, ~B), ~C) =>
                   broadcast((a,b,c)->(~f)((~g)(a,b), c), ~A, ~B, ~C) # TODO use syntactic dot fusion in the RHS -- needs some fix
                  ]

array_optimize(expr) = Fixpoint(Postwalk(Chain(array_rules)))(unwrap(expr))
