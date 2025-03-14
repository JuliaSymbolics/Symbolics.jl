using DifferentiationInterface, DifferentiationInterfaceTest
using Symbolics: Symbolics
using Test


test_differentiation(
    AutoSymbolics(), default_scenarios(; include_constantified=true); logging=LOGGING
);

test_differentiation(
    AutoSymbolics(),
    default_scenarios(; include_normal=false, include_cachified=true);
    excluded=[:jacobian],  # TODO: figure out why this fails
    logging=LOGGING,
);

test_differentiation(
    AutoSparse(AutoSymbolics()),
    sparse_scenarios(; band_sizes=0:-1);
    sparsity=true,
    logging=LOGGING,
);
