using Test
using Symbolics
using Symbolics: value
using Symbolics: VariableDefaultValue, VariableConnectType, VariableUnit
using Unitful

vals = [1,2,3,4]
@variables x=1 xs[1:4]=vals ys[1:5]=1

@test getmetadata(x, VariableDefaultValue) === 1
@test getmetadata.(xs, (VariableDefaultValue,)) == vals
@test getmetadata.(ys, (VariableDefaultValue,)) == ones(Int, 5)

struct Flow end
u = u"m^3/s"
@variables begin
    x = [1, 2], [connect=Flow,unit=u]
    y = 2
end

@test getmetadata(x, Symbolics.VariableDefaultValue) == [1, 2]
@test getmetadata(x, Symbolics.VariableConnectType) == Flow
@test getmetadata(x, Symbolics.VariableUnit) == u
@test getmetadata(y, Symbolics.VariableDefaultValue) === 2

@variables begin
    x, [connect=Flow,unit=u]
    y = 2, [connect=Flow]
end

@test !hasmetadata(x, Symbolics.VariableDefaultValue)
@test getmetadata(x, Symbolics.VariableConnectType) == Flow
@test getmetadata(x, Symbolics.VariableUnit) == u
@test getmetadata(y, Symbolics.VariableDefaultValue) === 2
@test getmetadata(y, Symbolics.VariableConnectType) == Flow

a = Symbolics.rename(value(x), :a)
@test !hasmetadata(x, Symbolics.VariableDefaultValue)
@test getmetadata(x, Symbolics.VariableConnectType) == Flow
@test getmetadata(x, Symbolics.VariableUnit) == u

@variables t x(t)=1 [connect=Flow,unit=u]

@test getmetadata(x, Symbolics.VariableDefaultValue) == 1
@test getmetadata(x, Symbolics.VariableConnectType) == Flow
@test getmetadata(x, Symbolics.VariableUnit) == u

a = Symbolics.rename(value(x), :a)
@test getmetadata(a, Symbolics.VariableDefaultValue) == 1
@test getmetadata(a, Symbolics.VariableConnectType) == Flow
@test getmetadata(a, Symbolics.VariableUnit) == u
