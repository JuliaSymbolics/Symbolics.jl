module SymbolicsHypergeometricFunctionsExt

using Symbolics
using HypergeometricFunctions

using Symbolics: @register_symbolic, @register_derivative

# Register the main hypergeometric functions as symbolic
# These are the fixed-arity versions that work well with Symbolics

# ₁F₁(a, b, z) - Kummer's confluent hypergeometric function
@register_symbolic HypergeometricFunctions._₁F₁(a, b, z)

# ₂F₁(a, b, c, z) - Gauss hypergeometric function
@register_symbolic HypergeometricFunctions._₂F₁(a, b, c, z) false

# ₃F₂(a₁, a₂, a₃, b₁, b₂, z) - Generalized hypergeometric function
@register_symbolic HypergeometricFunctions._₃F₂(a₁, a₂, a₃, b₁, b₂, z) false

# Register derivatives with respect to z (the last argument)
# These follow the standard differentiation formulas for hypergeometric functions:
#   d/dz ₁F₁(a, b, z) = (a/b) ₁F₁(a+1, b+1, z)
#   d/dz ₂F₁(a, b, c, z) = (ab/c) ₂F₁(a+1, b+1, c+1, z)
#   d/dz ₃F₂(a₁, a₂, a₃, b₁, b₂, z) = (a₁a₂a₃)/(b₁b₂) ₃F₂(a₁+1, a₂+1, a₃+1, b₁+1, b₂+1, z)

@register_derivative HypergeometricFunctions._₁F₁(a, b, z) 3 (a / b) * HypergeometricFunctions._₁F₁(a + 1, b + 1, z)

@register_derivative HypergeometricFunctions._₂F₁(a, b, c, z) 4 (a * b / c) * HypergeometricFunctions._₂F₁(a + 1, b + 1, c + 1, z)

@register_derivative HypergeometricFunctions._₃F₂(a₁, a₂, a₃, b₁, b₂, z) 6 (a₁ * a₂ * a₃) / (b₁ * b₂) * HypergeometricFunctions._₃F₂(a₁ + 1, a₂ + 1, a₃ + 1, b₁ + 1, b₂ + 1, z)

end
