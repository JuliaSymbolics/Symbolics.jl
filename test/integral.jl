using Symbolics, Test

@variables x y

I1 = Integral(x , y , 0)
I2 = Integral(x , y , 0)

@test I1 == I2
