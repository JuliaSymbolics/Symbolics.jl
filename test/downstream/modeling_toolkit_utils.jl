using Symbolics: diff2term
using ModelingToolkit
using ModelingToolkit: equivalent, get_unit, value, VariableUnit, UnitfulUnitCheck
using DynamicQuantities: @u_str as @du_str
using Unitful: @u_str as @uu_str
using Test

@variables t [unit = du"s"] k(t)
D = Differential(t)

@test equivalent(
    du"s^-1",
    get_unit(
        diff2term(
        value(D(k)),
        Dict(VariableUnit => get_unit(D(k)))
    )
    )
)

@test equivalent(
    du"s^-2",
    get_unit(
        diff2term(
        value(D(D(k))),
        Dict(VariableUnit => get_unit(D(D(k))))
    )
    )
)

@variables t [unit = uu"s"] k(t)
D = Differential(t)

@test equivalent(
    uu"s^-1",
    UnitfulUnitCheck.get_unit(
        diff2term(
        value(D(k)),
        Dict(VariableUnit => UnitfulUnitCheck.get_unit(D(k)))
    )
    )
)

@test equivalent(
    uu"s^-2",
    UnitfulUnitCheck.get_unit(
        diff2term(
        value(D(D(k))),
        Dict(VariableUnit => UnitfulUnitCheck.get_unit(D(D(k))))
    )
    )
)

@variables x y t
D = Differential(t)
Symbolics.diff2term(D(x))
Symbolics.diff2term(D(sqrt(sqrt(sqrt(+(x, y))))))
Symbolics.diff2term(D(x + y))
