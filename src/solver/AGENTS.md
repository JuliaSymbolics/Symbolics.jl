# Exact polynomial solving

- Keep elimination, exact division, and polynomial GCD coefficients in arbitrary-precision integer or rational arithmetic. Symbolic expression coefficient types alone do not guarantee that a dependency preserves precision internally.
- Return reduced rational functions. Validate specialization where an intermediate pivot vanishes but the full coefficient matrix remains nonsingular.
- Check parametric solver changes against an independent exact matrix solve and substitution into the original equations, including seeded 2×2 and 3×3 systems. Exercise the fast path directly as well as through `symbolic_solve`, so a fallback cannot conceal a broken optimization.
- Run all of `test/solver.jl` and `GROUP=Core` before publishing solver changes. Preserve the existing issue reproducer's performance assertions.
