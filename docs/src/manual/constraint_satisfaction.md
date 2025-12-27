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

[SymbolicSMT.jl](https://github.com/JuliaSymbolics/SymbolicSMT.jl) provides constraint satisfaction and theorem proving capabilities through the Z3 Theorem Prover. It allows resolving symbolic expressions under boolean constraints.

### Key Features

- **Z3 Integration**: Uses the powerful Z3 SMT (Satisfiability Modulo Theories) solver
- **Symbolic Constraints**: Work with symbolic mathematical expressions as constraints
- **Expression Resolution**: Resolve expressions to boolean values when provable under constraints
- **Boolean Logic**: Handle complex boolean formulas and logical relationships

### Basic Usage

```julia
using Symbolics, SymbolicSMT

@variables x::Integer

# Define constraints on a single variable
constraints = Constraints([x > 5])

# Resolve expressions under constraints
resolve(x > 0, constraints)   # true (provably true since x > 5 implies x > 0)
resolve(x < 0, constraints)   # false (provably false)
resolve(x > 10, constraints)  # x > 10 (cannot determine, returns original expression)

# Check if an expression is satisfiable under constraints
issatisfiable(x > 100, constraints)  # true - x could be 101

# Check if an expression is provably always true under constraints
isprovable(x > 0, constraints)  # true - always satisfied when x > 5
```

### Example Applications

#### Single Variable Constraints

```julia
using Symbolics, SymbolicSMT

@variables n::Integer

# Define bounds on a variable
constraints = Constraints([n >= 0, n <= 100])

# Check properties
isprovable(n >= 0, constraints)   # true (given as constraint)
issatisfiable(n > 50, constraints)  # true
issatisfiable(n > 200, constraints) # false (violates n <= 100)
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

### Example: Logic Puzzles

```julia
using Symbolics, SymbolicSMT

# Example: Solving a simple logic puzzle
# If A implies B, and B implies C, and A is true, then C must be true

@variables A::Bool B::Bool C::Bool

# Define implication: A ⟹ B is equivalent to (!A) | B
⟹(p, q) = (!p) | q

constraints = Constraints([
    A ⟹ B,  # A implies B
    B ⟹ C,  # B implies C
    A       # A is true
])

# Check if C is provably true under these constraints
isprovable(C, constraints)  # true
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

### Boolean Constraint Systems

SymbolicSMT.jl excels at boolean constraint systems with multiple variables:

```julia
using Symbolics, SymbolicSMT

@variables P::Bool Q::Bool R::Bool S::Bool

# Define implication helper
⟹(p, q) = (!p) | q

# Model a complex logical system
constraints = Constraints([
    P ⟹ Q,      # P implies Q
    Q ⟹ R,      # Q implies R
    R ⟹ S,      # R implies S
    P           # P is true
])

# Prove properties about the system
isprovable(S, constraints)      # true - S must be true
isprovable(Q & R, constraints)  # true - both Q and R must be true
issatisfiable(!Q, constraints)  # false - Q cannot be false
```

### Combining Boolean and Integer Constraints

```julia
using Symbolics, SymbolicSMT

@variables flag::Bool count::Integer

# Mixed constraints
constraints = Constraints([flag, count > 0])

# Check properties
isprovable(flag, constraints)      # true
isprovable(count > 0, constraints) # true
issatisfiable(count > 100, constraints) # true - count could be 101
```

!!! note "Current Limitations"
    SymbolicSMT.jl currently has limited support for multi-variable arithmetic expressions
    (e.g., `x + y > 0`). For best results, use boolean logic and single-variable
    arithmetic constraints. See the [SymbolicSMT.jl repository](https://github.com/JuliaSymbolics/SymbolicSMT.jl)
    for the latest updates on supported features.

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

!!! warning "Experimental Status"
    SymbolicSMT.jl is experimental software. Use it primarily for research and exploration. Always validate results through multiple approaches when possible.