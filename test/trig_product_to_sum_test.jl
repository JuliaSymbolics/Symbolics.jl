using Test
using Symbolics

@testset "Trigonometric Product-to-Sum Rules" begin
    @variables t ω φ ψ A B
    
    # Test basic cos(A) * cos(B) rule
    @test isequal(trigexpand(cos(A) * cos(B)), (cos(A - B) + cos(A + B)) / 2)
    
    # Test basic sin(A) * sin(B) rule  
    @test isequal(trigexpand(sin(A) * sin(B)), (cos(A - B) - cos(A + B)) / 2)
    
    # Test basic sin(A) * cos(B) rule
    @test isequal(trigexpand(sin(A) * cos(B)), (sin(A + B) + sin(A - B)) / 2)
    
    # Test basic cos(A) * sin(B) rule
    @test isequal(trigexpand(cos(A) * sin(B)), (sin(A + B) - sin(A - B)) / 2)
    
    # Test the specific case from issue #1644
    original = cos(ω*t + φ) * cos(ω*t + φ - ψ)
    expanded = trigexpand(original)
    expected = (cos(ψ) + cos(2*(ω*t + φ) - ψ)) / 2
    
    # Simplify both to compare
    expanded_simplified = simplify(expanded)
    expected_simplified = simplify(expected)
    
    @test isequal(expanded_simplified, expected_simplified)
    
    # Test that constants are preserved
    original_with_constants = 5 * cos(A) * cos(B)
    expanded_with_constants = trigexpand(original_with_constants)
    expected_with_constants = 5 * (cos(A - B) + cos(A + B)) / 2
    
    @test isequal(expanded_with_constants, expected_with_constants)
    
    println("All trigonometric product-to-sum tests passed!")
end