using LinearAlgebra, Symbolics, Test
@variables t, x(t), y(t), z(t)
@test Symbolics.islinear(x + y,[x,y])
@test Symbolics.islinear(x,[x,y])
@test Symbolics.islinear(y,[x,y])
@test !Symbolics.islinear(z,[x,y])

@test Symbolics.isaffine(x + y,[x,y])
@test Symbolics.isaffine(x,[x,y])
@test Symbolics.isaffine(y,[x,y])
@test Symbolics.isaffine(z,[x,y])
@test Symbolics.isaffine(x + y + z,[x,y])

@test Symbolics.isaffine(x + z * y,[x,y])
@test Symbolics.islinear(x + z * y,[x,y])
@test Symbolics.islinear(z * x + z * y,[x,y])
@test Symbolics.islinear(z * (x + y),[x,y])
@test Symbolics.isaffine(z * (x + y),[x,y])

@test Symbolics.isaffine(ifelse(x < 1, y, z), [z])
@test !Symbolics.isaffine(ifelse(x < 1, x, z), [x])

@test !Symbolics.isaffine(ifelse(y < 2, ifelse(x < 1, 0, 1), y), [x])

@variables w[1:2]
A = [1.0 2; 3 4]
@test Symbolics.islinear(sum(w), w)
@test Symbolics.islinear(sum(A * w), w)
@test Symbolics.islinear(sum(A * w .+ sum(w)), w)
@test Symbolics.islinear(sum(A * w), [w])
@test Symbolics.islinear(sum(w), [Symbolics.unwrap(w)])
@test Symbolics.islinear(sum(w) + z, [w, z])
@test Symbolics.islinear(sum(A * w) + z, [w, z])
@test Symbolics.islinear(sum(A * w), Symbolics.scalarize(w))
@test !Symbolics.islinear(sum(w) + 1, w)
@test !Symbolics.islinear(sum(w .^ 2), w)
@test !Symbolics.islinear(sum(w .^ 2), Symbolics.scalarize(w))
@test !Symbolics.islinear(prod(w), w)
@test !Symbolics.islinear(sum(w)^2, w)
@test !Symbolics.islinear(sum(sin.(w)), w)

v = sum(w)
@test Symbolics.islinear(2v, [v])
@test Symbolics.islinear(2(A * w)[1], [(A * w)[1]])
@test Symbolics.islinear(sum(A) * z, [z])

@variables C[1:2, 1:2]
@test !Symbolics.islinear(sum(C * w), [w; vec(C)])
@variables M[1:2, 1:2]
@test Symbolics.islinear(sum(M * w), w)
@test Symbolics.islinear(sum(M) * z, [z])
@test Symbolics.islinear(0, Num[])
@test !Symbolics.islinear(1, Num[])

@variables B[1:7, 1:7] p[1:5000] q
ex = det(B) * q
@test Symbolics.islinear(ex, [q])
Symbolics.islinear(ex, [q])
@test (@allocated Symbolics.islinear(ex, [q])) < 1_000_000
@test Symbolics.islinear(p[1] * q, [q])

function full_indexed_islinear_time(n)
    @variables scaling_x[1:n]
    vars = collect(Symbolics.scalarize(scaling_x))
    ex = vars[1] + vars[end]
    Symbolics.islinear(ex, vars)
    return minimum(@elapsed(Symbolics.islinear(ex, vars)) for _ in 1:3)
end

small_n_time = full_indexed_islinear_time(128)
large_n_time = full_indexed_islinear_time(512)
@test large_n_time < 10 * small_n_time
