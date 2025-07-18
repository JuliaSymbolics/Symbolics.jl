using Symbolics, Test
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
