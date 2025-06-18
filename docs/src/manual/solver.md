# Solver

The main symbolic solver for Symbolics.jl is `symbolic_solve`. Symbolic solving
means that it only uses symbolic (algebraic) methods and outputs exact solutions.

```@docs
Symbolics.symbolic_solve
```

One other symbolic solver is `symbolic_linear_solve` which is limited compared to 
`symbolic_solve` as it only solves linear equations.
```@docs
Symbolics.symbolic_linear_solve
```

`symbolic_solve` only supports symbolic, i.e. non-floating point computations, and thus prefers equations
where the coefficients are integer, rational, or symbolic. Floating point coefficients are transformed into
rational values and BigInt values are used internally with a potential performance loss, and thus it is recommended
that this functionality is only used with floating point values if necessary. In contrast, `symbolic_linear_solve`
directly handles floating point values using standard factorizations.

### More technical details and examples

#### Technical details

The `symbolic_solve` function uses 4 hidden solvers in order to solve the user's input. Its base,
`solve_univar`, uses analytic solutions up to polynomials of degree 4 and factoring as its method
for solving univariate polynomials. The function's `solve_multipoly` uses GCD on the input polynomials then throws passes the result
to `solve_univar`. The function's `solve_multivar` uses Groebner basis and a separating form in order to create linear equations in the
input variables and a single high degree equation in the separating variable [^1]. Each equation resulting from the basis is then passed
to `solve_univar`. We can see that essentially, `solve_univar` is the building block of `symbolic_solve`. If the input is not a valid polynomial and can not be solved by the algorithm above, `symbolic_solve` passes
it to `ia_solve`, which attempts solving by attraction and isolation [^2]. This only works when the input is a single expression
and the user wants the answer in terms of a single variable. Say `log(x) - a == 0` gives us `[e^a]`.

```@docs
Symbolics.solve_univar
Symbolics.solve_multivar
Symbolics.ia_solve
Symbolics.ia_conditions!
Symbolics.is_periodic
Symbolics.fundamental_period
```

#### Nice examples

```@example solver
using Symbolics, Nemo;
@variables x;
Symbolics.symbolic_solve(9^x + 3^x ~ 8, x)
```

```@example solver
@variables x y z;
Symbolics.symbolic_linear_solve(2//1*x + y - 2//1*z ~ 9//1*x, 1//1*x)
```

```@example solver
using Groebner;
@variables x y z;

eqs = [x^2 + y + z - 1, x + y^2 + z - 1, x + y + z^2 - 1]
Symbolics.symbolic_solve(eqs, [x,y,z])
```

### Feature completeness

- [x] Linear and polynomial equations
- [x] Systems of linear and polynomial equations
- [x] Some transcendental functions
- [x] Systems of linear equations with parameters (via `symbolic_linear_solve`)
- [ ] Equations with radicals
- [x] Systems of polynomial equations with parameters and positive dimensional systems
- [ ] Inequalities

### Expressions we can not solve (but aim to)
```
# Mathematica

In[1]:= Reduce[x^2 - x - 6 > 0, x]
Out[1]= x < -2 || x > 3

In[2]:= Reduce[x+a > 0, x]
Out[2]= a \[Element] Reals && x > -a

In[3]:= Solve[x^(x)  + 3 == 0, x]
Out[3]= {{x -> (I \[Pi] + Log[3])/ProductLog[I \[Pi] + Log[3]]}}
```

# References

[^1]: [Rouillier, F. Solving Zero-Dimensional Systems Through the Rational Univariate Representation. AAECC 9, 433–461 (1999).](https://doi.org/10.1007/s002000050114)
[^2]: [R. W. Hamming, Coding and Information Theory, ScienceDirect, 1980](https://www.sciencedirect.com/science/article/pii/S0747717189800070).
