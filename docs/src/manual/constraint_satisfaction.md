# Constraint Satisfaction

Symbolic constraint satisfaction is a powerful technique for solving problems where you need to find values for symbolic variables that satisfy a set of constraints. Symbolics.jl provides constraint satisfaction capabilities through integration with SAT solvers and theorem provers.

!!! note
    Constraint satisfaction in symbolic computation is an active area of research. The available tools provide experimental functionality that continues to evolve.

## What is Constraint Satisfaction?

A Constraint Satisfaction Problem (CSP) consists of:
- A set of variables
- A domain for each variable (the possible values it can take)
- A set of constraints that define relationships between variables

The goal is to find an assignment of values to variables that satisfies all constraints simultaneously.

In symbolic computation, constraint satisfaction extends beyond simple numerical domains to handle symbolic expressions, boolean formulas, and mathematical relationships.

## SymbolicSMT.jl

SymbolicSMT.jl extends SymbolicUtils expression simplification with theorem proving capabilities through the Z3 Theorem Prover. It allows adding boolean constraints to the symbolic simplification process.

### Key Features

- **Z3 Integration**: Uses the powerful Z3 SMT (Satisfiability Modulo Theories) solver
- **Symbolic Constraints**: Work with symbolic mathematical expressions as constraints
- **Satisfiability Checking**: Determine if expressions can be satisfied under constraints
- **Provability Checking**: Prove or disprove statements under constraints
- **Boolean Logic**: Handle complex boolean formulas and logical relationships

### Basic Usage

```@example constraint_sat
using Symbolics, SymbolicSMT

@variables x::Real y::Real

# Define constraints
constraints = Constraints([x > 0, y > 0])

# Check if an expression is satisfiable under constraints
issatisfiable(x + y > 1, constraints)
```

```@example constraint_sat
# Check if an expression is provable (always true) under constraints
isprovable(x >= 0, constraints)
```

```@example constraint_sat
# Check if x can be negative (it cannot, given x > 0)
issatisfiable(x < 0, constraints)
```

### Example Applications

#### Bound Checking

```@example constraint_sat
@variables a::Integer b::Integer

# Define bounds
bounds = Constraints([a >= 0, a <= 10, b >= 0, b <= 10])

# Can the sum exceed 15?
issatisfiable(a + b > 15, bounds)
```

```@example constraint_sat
# Is the sum always at most 20?
isprovable(a + b <= 20, bounds)
```

#### Quadratic Constraints

```@example constraint_sat
@variables p::Integer q::Integer

# Circle-like constraint
circle = Constraints([p^2 + q^2 < 25])

# Can p be 3?
issatisfiable(p == 3, circle)
```

```@example constraint_sat
# Can p be 5? (No, because 5^2 = 25 is not < 25)
issatisfiable(p == 5, circle)
```

## SAT Solving in Julia

### What is SAT?

Boolean Satisfiability Testing (SAT) is the problem of determining whether there exists an assignment of truth values to boolean variables that makes a given boolean formula evaluate to true. SAT is fundamental to many computational problems and is NP-complete.

### Applications of SAT Solvers

SAT solvers have wide applications in:

- **Software Verification**: Model checking, program analysis
- **Hardware Verification**: Circuit verification and testing
- **Planning and Scheduling**: Resource allocation, job scheduling
- **Artificial Intelligence**: Automated reasoning, knowledge representation
- **Cryptography**: Cryptanalysis and security analysis
- **Operations Research**: Optimization problems
- **Combinatorial Problems**: Sudoku, graph coloring, constraint puzzles

### SAT in Symbolic Computation

In symbolic computation, SAT solvers enable:

1. **Constraint Solving**: Finding symbolic solutions to constraint systems
2. **Theorem Proving**: Automatically proving or disproving mathematical statements
3. **Expression Simplification**: Simplifying expressions under logical constraints
4. **Satisfiability Checking**: Determining if constraint systems have solutions

### Boolean Constraints

```@example constraint_sat
@variables flag1::Bool flag2::Bool

# Boolean constraints
bool_constraints = Constraints([flag1, flag2])

# Both flags are true, so their AND is satisfiable
issatisfiable(flag1 & flag2, bool_constraints)
```

```@example constraint_sat
# Can flag1 be false? No, it's constrained to be true
issatisfiable(!flag1, bool_constraints)
```

## SMT Solvers

Satisfiability Modulo Theories (SMT) solvers extend SAT solvers to handle richer mathematical structures:

- **Arithmetic**: Integer and real number constraints
- **Arrays**: Array indexing and element constraints
- **Bit-vectors**: Fixed-width integer arithmetic
- **Strings**: String manipulation and pattern matching
- **Algebraic Data Types**: Custom mathematical structures

### Z3 Theorem Prover

Z3 is Microsoft Research's high-performance SMT solver that SymbolicSMT.jl uses:

- **Multiple Theories**: Supports arithmetic, bit-vectors, arrays, and more
- **Decision Procedures**: Efficient algorithms for theory-specific reasoning
- **Model Generation**: Can provide example solutions when constraints are satisfiable
- **Proof Generation**: Can generate proofs of unsatisfiability

## Advanced Constraint Satisfaction

### Multi-Variable Arithmetic

SymbolicSMT.jl supports multi-variable arithmetic expressions:

```@example constraint_sat
@variables m::Integer n::Integer

# Multi-variable constraints
multi_constraints = Constraints([m >= 1, n >= 1])

# Check multi-variable expressions - can sum be <= 0?
issatisfiable(m + n <= 0, multi_constraints)  # false - m + n >= 2 when m,n >= 1
```

```@example constraint_sat
# Is sum always >= 2?
isprovable(m + n >= 2, multi_constraints)  # true - always satisfied when m,n >= 1
```

### Power Operators

```@example constraint_sat
@variables t::Integer

power_constraints = Constraints([t >= 0, t <= 3])

# Can t^2 be at most 9?
issatisfiable(t^2 <= 9, power_constraints)
```

```@example constraint_sat
# Is t^2 always <= 9 when t is in [0,3]?
isprovable(t^2 <= 9, power_constraints)
```

### Expression Resolution

The `resolve` function attempts to determine if an expression is provably true or false:

```@example constraint_sat
@variables v::Integer

resolve_constraints = Constraints([v > 5])

# v > 0 is provably true when v > 5
resolve(v > 0, resolve_constraints)
```

```@example constraint_sat
# v < 0 is provably false when v > 5
resolve(v < 0, resolve_constraints)
```

```@example constraint_sat
# v > 10 cannot be determined - returns the expression
resolve(v > 10, resolve_constraints)
```

## Integration with Symbolic Computation

Constraint satisfaction integrates naturally with other symbolic computation features:

- **Symbolic Differentiation**: Optimize constraint systems using gradients
- **Symbolic Integration**: Handle constraints involving integrals
- **Expression Manipulation**: Simplify complex constraint expressions
- **Code Generation**: Generate efficient constraint checking code

## Further Reading

- [SymbolicSMT.jl Repository](https://github.com/JuliaSymbolics/SymbolicSMT.jl)
- [Z3 Theorem Prover](https://github.com/Z3Prover/z3)
- [Satisfiability.jl](https://github.com/dpsanders/SatisfiabilityInterface.jl) - Alternative SAT interface for Julia
- [SMT-LIB Standard](http://smtlib.cs.uiowa.edu/) - Standard format for SMT solvers

