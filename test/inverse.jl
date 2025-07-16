using Symbolics

@test inverse(sin) == left_inverse(sin) == right_inverse(sin) == asin
@test inverse(asin) == left_inverse(asin) == right_inverse(asin) == sin
@test has_inverse(sin) && has_left_inverse(sin) && has_right_inverse(sin)
fn = left_inverse(sqrt)
@test right_inverse(fn) == sqrt
@test_throws MethodError inverse(sqrt)
@test_throws MethodError right_inverse(sqrt)
@test !has_inverse(fn) && !has_left_inverse(fn)
@test !has_inverse(sqrt) && !has_right_inverse(sqrt)
@test has_inverse(sin ∘ cos)
@test !has_inverse(sin ∘ sqrt)
@test has_left_inverse(sin ∘ sqrt)
@test inverse(sin ∘ cos) == acos ∘ asin
@test inverse(inverse(sin ∘ cos)) == sin ∘ cos
@test right_inverse(left_inverse(sin ∘ sqrt)) == sin ∘ sqrt
@test inverse(sin ∘ cos ∘ tan) == atan ∘ (acos ∘ asin)
