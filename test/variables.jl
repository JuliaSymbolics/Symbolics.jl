using Test
using Symbolics
using Symbolics: DefaultValue

vals = [1,2,3,4]
@variables x=1 xs[1:4]=vals ys[1:5]=1

@test getmetadata(x, DefaultValue) === 1
@test getmetadata.(xs, (DefaultValue,)) == vals
@test getmetadata.(ys, (DefaultValue,)) == ones(Int, 5)
