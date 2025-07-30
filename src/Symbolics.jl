"""
$(DocStringExtensions.README)
"""
module Symbolics

using PrecompileTools

import PrecompileTools: @recompile_invalidations

@recompile_invalidations begin
    import CommonWorldInvalidations
end

using DocStringExtensions, Markdown

using LinearAlgebra

using Primes

using Reexport

using Setfield

import DomainSets: Domain, DomainSets

using TermInterface
import TermInterface: maketerm, iscall, operation, arguments, metadata

import SymbolicUtils: Term, Add, Mul, Pow, Sym, Div, BasicSymbolic,
FnType, @rule, Rewriters, substitute, symtype,
promote_symtype, isadd, ismul, ispow, isterm, issym, isdiv

using SymbolicUtils.Code

import SymbolicUtils.Rewriters: Chain, Prewalk, Postwalk, Fixpoint

import SymbolicUtils.Code: toexpr

import ArrayInterface
using RuntimeGeneratedFunctions
using SciMLBase
import MacroTools

using SymbolicIndexingInterface

import SymbolicLimits

using ADTypes: ADTypes

import OffsetArrays

@reexport using SymbolicUtils
RuntimeGeneratedFunctions.init(@__MODULE__)

# re-export

export simplify, substitute

export Num
import MacroTools: splitdef, combinedef, postwalk, striplines
include("wrapper-types.jl")

include("num.jl")

include("rewrite-helpers.jl")
include("complex.jl")

"""
    substitute(expr, s; fold=true)

Performs the substitution on `expr` according to rule(s) `s`.
If `fold=false`, expressions which can be evaluated won't be evaluated.
# Examples
```jldoctest
julia> @variables t x y z(t)
4-element Vector{Num}:
    t
    x
    y
 z(t)
julia> ex = x + y + sin(z)
(x + y) + sin(z(t))
julia> substitute(ex, Dict([x => z, sin(z) => z^2]))
(z(t) + y) + (z(t) ^ 2)
julia> substitute(sqrt(2x), Dict([x => 1]); fold=false)
sqrt(2)
```
"""
substitute

export Equation, ConstrainedEquation
include("equations.jl")

export Inequality, ≲, ≳
include("inequality.jl")

import Bijections, DynamicPolynomials
export tosymbol, terms, factors
include("utils.jl")

using ConstructionBase
include("arrays.jl")

export @register_symbolic, @register_array_symbolic
include("register.jl")

using SparseArrays
export @variables, Variable
include("variable.jl")

function slog end; function ssqrt end; function scbrt end
include("linearity.jl")

using DiffRules, SpecialFunctions, NaNMath


export Differential, expand_derivatives, is_derivative

include("diff.jl")

export SymbolicsSparsityDetector

include("adtypes.jl")

export Difference, DiscreteUpdate

include("difference.jl")

export infimum, supremum
include("domains.jl")

export Integral
include("integral.jl")

include("array-lib.jl")

using LogExpFunctions
include("logexpfunctions-lib.jl")

include("linear_algebra.jl")
export symbolic_linear_solve, solve_for

include("groebner_basis.jl")
export groebner_basis, is_groebner_basis

include("taylor.jl")
export series, taylor, taylor_coeff

import Libdl
include("build_function.jl")
export build_function

import Distributions
include("extra_functions.jl")

using Latexify
using LaTeXStrings
include("latexify_recipes.jl")

using RecipesBase
include("plot_recipes.jl")

include("semipoly.jl")


include("parsing.jl")
export parse_expr_to_symbolic

include("error_hints.jl")
include("struct.jl")
include("operators.jl")

include("limits.jl")
export limit

# Hacks to make wrappers "nicer"
const NumberTypes = Union{AbstractFloat,Integer,Complex{<:AbstractFloat},Complex{<:Integer}}
(::Type{T})(x::SymbolicUtils.Symbolic) where {T<:NumberTypes} = throw(ArgumentError("Cannot convert Sym to $T since Sym is symbolic and $T is concrete. Use `substitute` to replace the symbolic unwraps."))

"""
    _safe_rational_simplify(expr; kw...)

A wrapper around SymbolicUtils.simplify that validates the result when dealing with 
rational coefficients to work around known bugs in polynomial simplification.
See issue #1583.
"""
function _safe_rational_simplify(expr::Union{Num, Complex{Num}}; kw...)
    original = expr
    simplified = wrap(SymbolicUtils.simplify(unwrap(expr); kw...))
    
    # Check if the expression involves rational coefficients
    if _has_rational_coeffs(unwrap(original))
        # Validate the simplification by testing with sample values
        if !_validate_simplification(original, simplified)
            # If validation fails, return the original expression
            return original
        end
    end
    
    return simplified
end

"""
    _has_rational_coeffs(expr)

Check if an expression contains rational number coefficients (///).
"""
function _has_rational_coeffs(expr)
    if isa(expr, Rational)
        return true
    elseif iscall(expr)
        return any(_has_rational_coeffs, arguments(expr))
    else
        return false
    end
end

"""
    _validate_simplification(original, simplified)

Validate that a simplification is mathematically correct by testing with random values.
"""
function _validate_simplification(original::Union{Num, Complex{Num}}, simplified::Union{Num, Complex{Num}})
    try
        # Get all variables in both expressions
        orig_vars = get_variables(original)
        simp_vars = get_variables(simplified)
        all_vars = unique([orig_vars; simp_vars])
        
        # If no variables, compare directly
        if isempty(all_vars)
            return isequal(original, simplified)
        end
        
        # Test with a few random values to check mathematical equivalence
        for _ in 1:3
            test_vals = Dict(var => rand(1:10) for var in all_vars)
            
            # Evaluate both expressions
            orig_val = substitute(original, test_vals)
            simp_val = substitute(simplified, test_vals)
            
            # Check if values are approximately equal (allowing for floating point errors)
            if !isapprox(Float64(orig_val), Float64(simp_val), rtol=1e-10)
                return false
            end
        end
        
        return true
    catch
        # If validation fails for any reason, be conservative and reject the simplification
        return false
    end
end

for T in [Num, Complex{Num}]
    @eval begin
        #(::Type{S})(x::$T) where {S<:Union{NumberTypes,AbstractArray}} = S(Symbolics.unwrap(x))::S

        SymbolicUtils.simplify(n::$T; kw...) = _safe_rational_simplify(n; kw...)
        SymbolicUtils.simplify_fractions(n::$T; kw...) = wrap(SymbolicUtils.simplify_fractions(unwrap(n); kw...))
        SymbolicUtils.expand(n::$T) = wrap(SymbolicUtils.expand(unwrap(n)))
        substitute(expr::$T, s::Pair; kw...) = wrap(substituter(s)(unwrap(expr); kw...)) # backward compat
        substitute(expr::$T, s::Vector; kw...) = wrap(substituter(s)(unwrap(expr); kw...))
        substitute(expr::$T, s::Dict; kw...) = wrap(substituter(s)(unwrap(expr); kw...))

        SymbolicUtils.Code.toexpr(x::$T) = SymbolicUtils.Code.toexpr(unwrap(x))

        SymbolicUtils.setmetadata(x::$T, t, v) = wrap(SymbolicUtils.setmetadata(unwrap(x), t, v))
        SymbolicUtils.getmetadata(x::$T, t) = SymbolicUtils.getmetadata(unwrap(x), t)
        SymbolicUtils.hasmetadata(x::$T, t) = SymbolicUtils.hasmetadata(unwrap(x), t)

        Broadcast.broadcastable(x::$T) = x
    end
    for S in [:(Symbolic{<:FnType}), :CallWithMetadata]
        @eval (f::$S)(x::$T, y...) = wrap(f(unwrap(x), unwrap.(y)...))
    end
end

for sType in [Pair, Vector, Dict]
    @eval substitute(expr::Arr, s::$sType; kw...) = wrap(substituter(s)(unwrap(expr); kw...))
end

# Symbolic solver
include("solver/preprocess.jl")
include("solver/nemo_stuff.jl")
include("solver/solve_helpers.jl")
include("solver/postprocess.jl")
include("solver/univar.jl")
include("solver/ia_helpers.jl")
include("solver/polynomialization.jl")
include("solver/attract.jl")
include("solver/ia_main.jl")
include("solver/main.jl")
include("solver/special_cases.jl")
export symbolic_solve

# Sympy Functions

"""
    symbolics_to_sympy(expr)

Converts a Symbolics.jl expression to a SymPy expression. It is represented in
[SymPy.jl](https://juliapy.github.io/SymPy.jl/dev/) format and thus can use all
of the wrapped functionality from SymPy.jl.

For conversion back to Symbolics, see `sympy_to_symbolics`.

# Arguments

- `expr`: A Symbolics.jl expression

# Example

```julia
using Symbolics
@variables x y
expr = x^2 + y
sympy_expr = symbolics_to_sympy(expr)
```
"""
function symbolics_to_sympy end

"""
    sympy_to_symbolics(sympy_expr, vars)
Converts a SymPy expression to Symbolics.jl.
# Arguments
- `sympy_expr`: SymPy expression.
- `vars`: List or dictionary of Symbolics variables.
# Example
```julia
@variables x y
sympy_expr = SymPy.Sym("x")^2 + SymPy.Sym("y")
symbolics_expr = sympy_to_symbolics(sympy_expr, [x, y])
```
"""
function sympy_to_symbolics end

"""
    sympy_linear_solve(A, b)
Solves linear system Ax = b using SymPy.
# Arguments
- `A`: Matrix of Symbolics expressions.
- `b`: Vector of Symbolics expressions.
# Returns
Vector of Symbolics solutions.
# Example
```julia
@variables x y
A = [1 2; 3 4]
b = [x, y]
sol = sympy_linear_solve(A, b)
```
"""
function sympy_linear_solve end

"""
    sympy_algebraic_solve(expr, var)
Solves algebraic equation(s) expr = 0 for var(s) using SymPy.
# Arguments
- `expr`: Symbolics expression or vector of expressions for a system of equations (linear or nonlinear).
- `var`: Symbolics variable or vector of variables to solve for.
# Returns
- For a single equation: List of Symbolics solutions.
- For a system: List of dictionaries mapping variables to solutions.
# Example
```julia
@variables x y
# Single equation
expr = x^2 - 4
sol = sympy_algebraic_solve(expr, x)  # Returns [2, -2]
# Nonlinear system
eqs = [x^2 + y^2 - 4, x - y]  # Circle and line
sol = sympy_algebraic_solve(eqs, [x, y])  # Returns [{x=>1, y=>1}, {x=>-1, y=>-1}]
```
"""
function sympy_algebraic_solve end

"""
    sympy_integrate(expr, var)
Computes indefinite integral of expr w.r.t. var using SymPy.
# Arguments
- `expr`: Symbolics expression.
- `var`: Symbolics variable.
# Returns
Symbolics integral.
# Example
```julia
@variables x
expr = x^2
result = sympy_integrate(expr, x)
```
"""
function sympy_integrate end

"""
    sympy_limit(expr, var, val)
Computes limit of expr as var approaches val using SymPy.
# Arguments
- `expr`: Symbolics expression.
- `var`: Symbolics variable.
- `val`: Symbolics expression or number.
# Returns
Symbolics limit.
# Example
```julia
@variables x
expr = 1/x
result = sympy_limit(expr, x, 0)
```
"""
function sympy_limit end

"""
    sympy_simplify(expr)
Simplifies a Symbolics expression using SymPy.
# Arguments
- `expr`: Symbolics expression.
# Returns
Simplified Symbolics expression.
# Example
```julia
@variables x
expr = x^2 + 2x^2
result = sympy_simplify(expr)
```
"""
function sympy_simplify end

"""
    sympy_ode_solve(expr, func, var)
Solves ODE expr = 0 for function func w.r.t. var using SymPy.
# Arguments
- `expr`: Symbolics expression representing ODE (set to 0).
- `func`: Symbolics function (e.g., f(x)).
- `var`: Independent Symbolics variable.
# Returns
Symbolics solution(s).
# Example
```julia
@variables x
@syms f(x)
expr = Symbolics.Derivative(f, x) - 2*f
sol = sympy_ode_solve(expr, f, x)  # Returns C1*exp(2*x)
```
"""
function sympy_ode_solve end

export symbolics_to_sympy, sympy_to_symbolics
export sympy_linear_solve, sympy_algebraic_solve, sympy_integrate, sympy_limit, sympy_simplify,sympy_ode_solve

function __init__()
    Base.Experimental.register_error_hint(TypeError) do io, exc
        if exc.expected == Bool && exc.got isa Num
            println(io,
                "\nA symbolic expression appeared in a Boolean context. This error arises in situations where Julia expects a Bool, like ")
            printstyled(io, "if boolean_condition", color = :blue)
            printstyled(
                io, "\t\t use ifelse(boolean_condition, then branch, else branch)\n",
                color = :green)
            printstyled(io, "x && y", color = :blue)
            printstyled(io, "\t\t\t\t use x & y\n", color = :green)
            printstyled(io, "boolean_condition ? a : b", color = :blue)
            printstyled(io, "\t use ifelse(boolean_condition, a, b)\n", color = :green)
            print(io,
                "but a symbolic expression appeared instead of a Bool. For help regarding control flow with symbolic variables, see https://docs.sciml.ai/ModelingToolkit/dev/basics/FAQ/#How-do-I-handle-if-statements-in-my-symbolic-forms?")
        end
    end
end

export inverse, left_inverse, right_inverse, @register_inverse, has_inverse, has_left_inverse, has_right_inverse
include("inverse.jl")

export rootfunction, left_continuous_function, right_continuous_function, @register_discontinuity
include("discontinuities.jl")

end # module
