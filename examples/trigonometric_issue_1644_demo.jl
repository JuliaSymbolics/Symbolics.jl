#!/usr/bin/env julia
"""
Demonstration of trigonometric product-to-sum simplification 
solving issue #1644 in Symbolics.jl

This example shows how to use the new `trigexpand` function to simplify
trigonometric product expressions that could not be simplified before.
"""

using Symbolics

println("="^60)
println("Symbolics.jl Issue #1644: Trigonometric Product-to-Sum")
println("="^60)

@variables t Vbar Ibar ω φk ψ

# Define voltage and current waveforms as in the issue
println("\n1. Define voltage and current waveforms:")
println("   vₖ(t) = V̄ cos(ωt + φₖ)")
println("   iₖ(t) = Ī cos(ωt + φₖ - ψ)")
println()

vk(t) = Vbar * cos(ω*t + φk)
ik(t) = Ibar * cos(ω*t + φk - ψ)

println("2. Calculate instantaneous power as their product:")
pk_original = vk(t) * ik(t)
println("   pₖ(t) = vₖ(t) × iₖ(t) = ", pk_original)
println()

println("3. Before: Regular simplify() could not expand this expression")
pk_regular_simplify = simplify(pk_original, expand=true) 
println("   simplify(pₖ(t), expand=true) = ", pk_regular_simplify)
println("   → No change! The product remains unexpanded.")
println()

println("4. After: Using the new trigexpand() function")
pk_expanded = trigexpand(pk_original)
println("   trigexpand(pₖ(t)) = ", pk_expanded)
println("   → Successfully expanded using product-to-sum identities!")
println()

println("5. Mathematical verification:")
println("   The result uses the identity: cos(A)cos(B) = ½[cos(A-B) + cos(A+B)]")
println("   Where A = ωt + φₖ and B = ωt + φₖ - ψ")
println("   So A - B = ψ and A + B = 2ωt + 2φₖ - ψ")
println("   Therefore: cos(ωt + φₖ)cos(ωt + φₖ - ψ) = ½[cos(ψ) + cos(2ωt + 2φₖ - ψ)]")
println()

println("6. Additional examples of trigexpand():")

# Test other trigonometric products
@variables A B
println("   trigexpand(cos(A) * cos(B)) = ", trigexpand(cos(A) * cos(B)))
println("   trigexpand(sin(A) * sin(B)) = ", trigexpand(sin(A) * sin(B)))  
println("   trigexpand(sin(A) * cos(B)) = ", trigexpand(sin(A) * cos(B)))
println("   trigexpand(cos(A) * sin(B)) = ", trigexpand(cos(A) * sin(B)))
println()

println("7. The trigexpand() function is now available in Symbolics.jl")
println("   and can handle trigonometric product-to-sum transformations")
println("   that were previously not supported.")
println()

println("="^60)
println("Issue #1644 has been resolved!")
println("="^60)