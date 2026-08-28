# Every entry is a call that once raised a dispatch `MethodError` (mostly "is ambiguous")
# on user-facing mixed symbolic/numeric operations. The exact result is not checked here;
# the point is that each call resolves to a method and returns.
using Symbolics, LinearAlgebra, SparseArrays, Test
using StaticArraysCore: SArray
using StaticArrays: SHermitianCompact
using Symbolics: SpecialFunctions

@variables x y z t
@variables A[1:2, 1:2] b[1:2] c[1:2] R[1:1, 1:2]
@variables cb[1:2]::Complex
D = Differential(t)
M = [1.0 2.0; 3.0 4.0]
Mn = collect(A)
v = [1.0, 2.0]
S = SArray{Tuple{2, 2}, Float64, 2, 4}((1.0, 2.0, 3.0, 4.0))
sv = SArray{Tuple{2}}((1.0, 2.0))
Dg = Diagonal([1.0, 2.0])
Dgn = Diagonal([x, y])
UT = UpperTriangular(M)
UTn = UpperTriangular([x y; 0 z])
sp = sparse(M)

const CALLS = [
    "Num(1 + 0im)" => () -> Num(1 + 0im),
    "promote_type(ComplexF64, Num)" => () -> promote_type(ComplexF64, Num),
    "[1 + 2im, x]" => () -> [1 + 2im, x],
    "(1 + 2im) == x" => () -> (1 + 2im) == x,
    "x == pi" => () -> x == pi,
    "isequal(x, 1 + 2im)" => () -> isequal(x, 1 + 2im),
    "(x + im * y)^true" => () -> (x + im * y)^true,
    "1 ~ 3" => () -> 1 ~ 3,
    "[x, y] ~ [1, 2]" => () -> [x, y] ~ [1, 2],
    "Dg \\ b" => () -> Dg \ b,
    "Dgn \\ b" => () -> Dgn \ b,
    "Diagonal(sv) \\ b" => () -> Diagonal(sv) \ b,
    "UT \\ b" => () -> UT \ b,
    "UTn \\ b" => () -> UTn \ b,
    "sp \\ b" => () -> sp \ b,
    "sp' \\ b" => () -> sp' \ b,
    "Symmetric(M) \\ b" => () -> Symmetric(M) \ b,
    "Bidiagonal \\ b" => () -> Bidiagonal([1.0, 2.0], [3.0], :U) \ b,
    "SymTridiagonal \\ b" => () -> SymTridiagonal([1.0, 2.0], [3.0]) \ b,
    "S \\ b" => () -> S \ b,
    "M \\ b" => () -> M \ b,
    "A \\ b" => () -> A \ b,
    "A \\ v" => () -> A \ v,
    "x \\ Dg" => () -> x \ Dg,
    "x \\ sv" => () -> x \ sv,
    "x \\ UT" => () -> x \ UT,
    "Dg / x" => () -> Dg / x,
    "Dgn / x" => () -> Dgn / x,
    "UT / x" => () -> UT / x,
    "Symmetric(M) / x" => () -> Symmetric(M) / x,
    "sparsevec(v) / x" => () -> sparsevec(v) / x,
    "view(sp, :, 1) / x" => () -> view(sp, :, 1) / x,
    "sv / x" => () -> sv / x,
    "S / x" => () -> S / x,
    "SHermitianCompact / x" => () -> SHermitianCompact{2}(SArray{Tuple{3}}((1.0, 2.0, 3.0))) / x,
    "x \\ SHermitianCompact" => () -> x \ SHermitianCompact{2}(SArray{Tuple{3}}((1.0, 2.0, 3.0))),
    "(0.0:0.5:1.0) / x" => () -> (0.0:0.5:1.0) / x,
    "BitVector / x" => () -> BitVector([1, 0]) / x,
    "A / Dgn" => () -> A / Dgn,
    "A / UTn" => () -> A / UTn,
    "A / M" => () -> A / M,
    "M / A" => () -> M / A,
    "collect(b)' / A" => () -> collect(b)' / A,
    "2 * A * b" => () -> 2 * A * b,
    "x * A * b" => () -> x * A * b,
    "A * b * 2" => () -> A * b * 2,
    "2 * A * A * b" => () -> 2 * A * A * b,
    "M * A * b" => () -> M * A * b,
    "A * M * b" => () -> A * M * b,
    "A * A * A * A * A" => () -> A * A * A * A * A,
    "S * b" => () -> S * b,
    "S * A" => () -> S * A,
    "Dgn * b" => () -> Dgn * b,
    "A * Dg" => () -> A * Dg,
    "UT * b" => () -> UT * b,
    "UTn * A" => () -> UTn * A,
    "A * UT" => () -> A * UT,
    "b * v'" => () -> b * v',
    "b * adjoint(reshape(v, 2, 1))" => () -> b * adjoint(reshape(v, 2, 1)),
    "v * R" => () -> v * R,
    "b' * A * b" => () -> b' * A * b,
    "A + S" => () -> A + S,
    "S + A" => () -> S + A,
    "map(+, b, b)" => () -> map(+, b, b),
    "map(+, b, b, b)" => () -> map(+, b, b, b),
    "map(+, b, sv)" => () -> map(+, b, sv),
    "map(+, sv, b)" => () -> map(+, sv, b),
    "map(+, b, v)" => () -> map(+, b, v),
    "map(+, v, b)" => () -> map(+, v, b),
    "mapreduce(*, +, b, b)" => () -> mapreduce(*, +, b, b),
    "mapreduce(*, +, b, v)" => () -> mapreduce(*, +, b, v),
    "mapreduce(*, +, sv, b)" => () -> mapreduce(*, +, sv, b),
    "x in 1:3" => () -> x in 1:3,
    "x in 1.0:0.5:3.0" => () -> x in 1.0:0.5:3.0,
    "x in sv" => () -> x in sv,
    "x in [x, y]" => () -> x in [x, y],
    "x in view([x, y], 1:2)" => () -> x in view([x, y], 1:2),
    "ifelse(x > 0, b, cb)" => () -> ifelse(x > 0, b, cb),
    "ldiv!(UTn, sparsevec)" => () -> ldiv!(UpperTriangular([x y; 0 z]), sparsevec(Num[x, y])),
    "ldiv!(UTn, view(v, 1:2))" => () -> ldiv!(UpperTriangular([x y; 0 z]), view(Num[x, y], 1:2)),
    "D * D" => () -> D * D,
    "D * identity" => () -> D * identity,
    "copysign(1, x)" => () -> copysign(1, x),
    "copysign(1.0, x)" => () -> copysign(1.0, x),
    "copysign(1 // 2, x)" => () -> copysign(1 // 2, x),
    "besselj(x, 1.0)" => () -> SpecialFunctions.besselj(x, 1.0),
    "polygamma(2, x)" => () -> SpecialFunctions.polygamma(2, x),
    "x^M" => () -> x^M,
    "M^x" => () -> M^x,
    "Mn^x" => () -> Mn^x,
    "Dg^x" => () -> Dg^x,
    "A^x" => () -> A^x,
    "T3[i, 1, j]" => () -> begin
        @variables T3[1:2, 1:2, 1:2] i::Int j::Int
        T3[i, 1, j]
    end,
]

@testset "dispatch smoke: $name" for (name, f) in CALLS
    @test (f(); true)
end
