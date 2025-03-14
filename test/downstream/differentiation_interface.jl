using DifferentiationInterface, DifferentiationInterfaceTest
using Symbolics: Symbolics
using Test

LOGGING = get(ENV, "CI", "false") == "false"

for backend in [AutoSymbolics(), AutoSparse(AutoSymbolics())]
    @test check_available(backend)
    @test check_inplace(backend)
end

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
