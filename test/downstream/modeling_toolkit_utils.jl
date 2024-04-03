using Symbolics: diff2term
using ModelingToolkit
using ModelingToolkit: equivalent, get_unit, value, VariableUnit
using Unitful: @u_str
using Test

@variables t [unit = u"s"] k(t)
D = Differential(t)

@test equivalent(
    u"1/s",
    get_unit(
        diff2term(
            value(D(k)),
            Dict(VariableUnit => get_unit(D(k)))
        )
    )
)

@test equivalent(
    u"1/s^2",
    get_unit(
        diff2term(
            value(D(D(k))),
            Dict(VariableUnit => get_unit(D(D(k))))
        )
    )
)
