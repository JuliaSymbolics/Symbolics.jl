# Solver
The main `symbolic` solver for Symbolics.jl is `symbolic_solve`. Symbolic solving
means that it only uses analytical methods and manually solves equations for input variables.
It uses no numerical methods and only outputs exact solutions.
```@docs
Symbolics.symbolic_solve
```
### More technical details
This function uses 4 solvers in order to solve the user's input. Its base,
`solve_univar`, uses analytic solutions up to polynomials of degree 4 and factoring as its method
for solving univariate polynomials. The function's `solve_multipoly` uses GCD on the input polynomials then throws passes the result
to `solve_univar`. The function's `solve_multivar` uses Groebner basis and a separating form in order to create linear equations in the
input variables and a single high degree equation in the separating variable. Each equation resulting from the basis is then passed
to `solve_univar`. We can see that essentially, `solve_univar` is the building block of `symbolic_solve`.if the input is not a valid polynomial and can not be solved by the algorithm above, `symbolic_solve` passes
it to `ia_solve`, which attempts solving by attraction and isolation [^1]. This only works when the input is a single expression
and the user wants the answer in terms of a single variable. Say `log(x) - a == 0` gives us `[e^a]`.

# References
[^1]: [R. W. Hamming, Coding and Information Theory, ScienceDirect, 1980](https://www.sciencedirect.com/science/article/pii/S0747717189800070).

