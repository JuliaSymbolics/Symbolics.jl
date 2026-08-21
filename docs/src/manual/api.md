# Public API

This page collects the public names that are part of the Symbolics.jl interface.
The linked docstrings describe the arguments, keyword arguments, return values, and
usage contracts for each API.

```@docs
Symbolics.@register_discontinuity
Symbolics.Inequality
Symbolics.Num
Symbolics._toexpr_metadata
Symbolics._toexpr_op
Symbolics.approximation_function
Symbolics.gather_factor
Symbolics.get_differential_vars
Symbolics.arguments
Symbolics.hasnode
Symbolics.filterchildren
Symbolics.replacenode
Symbolics.is_groebner_basis
Symbolics.left_continuous_function
Symbolics.majorization_function
Symbolics.minorization_function
Symbolics.partial_frac_decomposition
Symbolics.polynomial_coeffs
Symbolics.right_continuous_function
Symbolics.rootfunction
Symbolics.semilinear_form
Symbolics.semipolynomial_form
Symbolics.semiquadratic_form
Symbolics.substitute_in_deriv
Symbolics.substitute_in_deriv_and_depvar
Symbolics.symbolics_to_sympy_pythoncall
Symbolics.sympy_pythoncall_algebraic_solve
Symbolics.sympy_pythoncall_integrate
Symbolics.sympy_pythoncall_limit
Symbolics.sympy_pythoncall_linear_solve
Symbolics.sympy_pythoncall_simplify
Symbolics.sympy_pythoncall_to_symbolics
Symbolics.value
Symbolics.:≲
Symbolics.:≳
Symbolics.wrap
Symbolics.Symbolics
Symbolics.parse_vars
Symbolics.FixpointSubstituter
```

## Reexported dependency APIs

These names are implemented and documented by their owning packages and reexported
by Symbolics. The owner documentation is canonical.

```@docs
SymbolicUtils.substitute
```

The reexported `simplify`, `unwrap`, `ifelse_eager`, and `ifelse_branching` APIs
are documented by [SymbolicUtils](https://symbolicutils.juliasymbolics.org/api/).
