using Symbolics
using LinearAlgebra
using SparseArrays: sparse
using Test


@variables a,b,c,d,e,f,g,h,i

# test hashing
aa = a; # old a

@variables a

@test isequal(a, aa)
@test hash(a) == hash(aa)

# test some matrix operations don't throw errors
X = [0 b c; d e f; g h i]
@test iszero(simplify(det(X) - ((d * ((b * i) - (c * h))) + (g * ((b * f) - (c * e))))))
F = lu(X)
@test F.p == [2, 1, 3]
R = simplify.(F.L * F.U - X[F.p, :], polynorm=true)
@test iszero(R)
@test simplify.(F \ X) == I
@test Symbolics._solve(X, X) == I
inv(X)
qr(X)

X2 = [0 b c; 0 0 0; 0 h 0]
@test_throws SingularException lu(X2)
F2 = lu(X2, check=false)
@test F2.info == 1

# test operations with sparse arrays and Operations
# note `isequal` instead of `==` because `==` would give another Operation

# test that we can create a sparse array of Operation
Oarray = zeros(Num, 2,2)
Oarray[2,2] = a
@test isequal(sparse(Oarray), sparse([2], [2], [a]))

# test Operation * sparse
@test isequal(a * sparse([2], [2], [1]), sparse([2], [2], [a * 1]))

# test sparse{Operation} + sparse{Operation}
A = sparse([2], [2], [a])
B = sparse([2], [2], [b])
@test isequal(A + B, sparse([2], [2], [a+b]))

# test sparse{Operation} * sparse{Operation}
C = sparse([1, 2], [2, 1], [c, c])
D = sparse([1, 2], [2, 1], [d, d])

@test isequal(C * D, sparse([1,2], [1,2], [c * d, c * d]))

@variables t σ ρ β
@variables x(t) y(t) z(t)

Dx = Differential(x)
Dy = Differential(y)
Dz = Differential(z)

eqs = [σ*(y-x),
       x*(ρ-z)-y,
       x*y - β*z]
J = Num[Dx(eqs[1]) Dy(eqs[1]) Dz(eqs[1])
        Dx(eqs[2]) Dy(eqs[2]) Dz(eqs[2])
        Dx(eqs[3]) Dy(eqs[3]) Dz(eqs[3])]

J = expand_derivatives.(J)
using LinearAlgebra
luJ = lu(J,Val(false))

@variables M[1:2,1:2]
inv(M)

@variables b[1:2]
M = [1 0; 0 2]
M \ b
M \ reshape(b,2,1)
M = [1 1; 0 2]
M \ reshape(b,2,1)


M = [1 a; 0 2]
M \ b
M \ [1, 2]
