# Trigonometric Product-to-Sum Simplification Rules
# These rules implement the product-to-sum trigonometric identities

using SymbolicUtils: @rule, @acrule
using SymbolicUtils.Rewriters: Chain, Prewalk

"""
Additional trigonometric product-to-sum simplification rules to complement 
the existing SymbolicUtils trigonometric simplification.

These rules implement:
- cos(A) * cos(B) = (1/2) * [cos(A-B) + cos(A+B)]
- sin(A) * sin(B) = (1/2) * [cos(A-B) - cos(A+B)]  
- sin(A) * cos(B) = (1/2) * [sin(A+B) + sin(A-B)]
- cos(A) * sin(B) = (1/2) * [sin(A+B) - sin(A-B)]
"""
const TRIG_PRODUCT_TO_SUM_RULES = [
    # cos(A) * cos(B) = (1/2) * [cos(A-B) + cos(A+B)]
    @acrule(cos(~A) * cos(~B) => (cos(~A - ~B) + cos(~A + ~B)) / 2),
    
    # sin(A) * sin(B) = (1/2) * [cos(A-B) - cos(A+B)]
    @acrule(sin(~A) * sin(~B) => (cos(~A - ~B) - cos(~A + ~B)) / 2),
    
    # sin(A) * cos(B) = (1/2) * [sin(A+B) + sin(A-B)]
    @acrule(sin(~A) * cos(~B) => (sin(~A + ~B) + sin(~A - ~B)) / 2),
    
    # cos(A) * sin(B) = (1/2) * [sin(A+B) - sin(A-B)]
    @acrule(cos(~A) * sin(~B) => (sin(~A + ~B) - sin(~A - ~B)) / 2)
]

"""
    trigexpand(expr)

Apply trigonometric product-to-sum expansion rules to the expression.
This converts products of trigonometric functions to sums.

# Examples
```julia
@variables t ω φ ψ
expr = cos(ω*t + φ) * cos(ω*t + φ - ψ)
trigexpand(expr)  # Returns expanded form using product-to-sum identities
```
"""
function trigexpand(expr)
    # Apply the product-to-sum rules using SymbolicUtils rewriter
    rewriter = Prewalk(Chain(TRIG_PRODUCT_TO_SUM_RULES))
    result = rewriter(unwrap(expr))
    
    # Apply standard simplification after expansion
    simplified = SymbolicUtils.simplify(result)
    return wrap(simplified)
end

# Export the function for external use
export trigexpand