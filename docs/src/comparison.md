# Comparison of Julia's Symbolics.jl vs SymPy for Symbolic Computation

Symbolics.jl is a symbolic modeling language for Julia, built in
Julia. Its goal is very different from Sympy: it was made to support
symbolic-numerics, the combination of symbolic computing with numerical
methods to allow for extreme performance computing that would not be
possible without modifying the model. Because of this, Symbolics.jl
excels in many areas due to purposeful design decisions:

- Performance: Symbolics.jl is built in Julia, whereas SymPy was
  built in Python. Thus, the performance bar for Symbolics.jl is
  much higher. Symbolics.jl started because SymPy was far too
  slow and SymEngine was far too inflexible for the projects they were
  doing. Performance is key to Symbolics.jl. If you find any
  performance issues, please file an issue.
- `build_function`: `lambdify` is “fine” for some people, but if you're building
  a super fast MPI-enabled Julia/C/Fortran simulation code, having a
  function that hits the Python interpreter is less than optimal. By
  default, `build_function` builds fast JIT-compiled functions due
  to being in Julia. However, it has support for things like static
  arrays, non-allocating functions via mutation, fast functions on
  sparse matrices and arrays of arrays, etc.: all core details of
  doing high performance computing.
- Parallelism: Symbolics.jl has pervasive parallelism. The
  symbolic simplification via [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl)
  has built-in parallelism, Symbolics.jl builds functions that
  parallelize across threads. Symbolics.jl is compatible with GPU libraries like CUDA.jl.
- Extensible: Symbolics.jl and its underlying tools are written in
  pure Julia. Want to add new or better simplification rules? Add some Julia code!
  Need to add new derivatives? Add some Julia code! You get the picture. Breaking
  down these barriers makes it easier for the user to tailor the program to their
  needs and accelerates the development of the library.
- Deep integration with the Julia ecosystem: Symbolics.jl's integration
  with neural networks is not the only thing that's deep. Symbolics.jl
  is built with the same philosophy as other SciML packages, eschewing
  “monorepos” for a distributed development approach that ties together
  the work of many developers. The differentiation parts utilize tools
  from automatic differentiation libraries, all linear algebra functionality
  comes from tracing Julia Base itself, symbolic rewriting (simplification
  and substitution) comes from
  [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl),
  parallelism comes from Julia Base libraries, Dagger.jl, etc.
  SciML Tools like
  [DataDrivenDiffEq.jl](https://datadriven.sciml.ai/dev/) can reconstruct
  symbolic expressions from neural networks and data, while
  [NeuralPDE.jl](https://github.com/SciML/NeuralPDE.jl)
  can automatically solve partial differential equations from symbolic
  descriptions using physics-informed neural networks.
  The list keeps going. All told, by design Symbolics.jl's development
  moves fast because it's effectively using the work of hundreds of
  Julia developers, allowing it to grow fast.
- While Symbolics.jl has many features missing from SymPy, it does not superset
  SymPy's functionality. For a list of missing features, see [this issue](https://github.com/JuliaSymbolics/Symbolics.jl/issues/59).
