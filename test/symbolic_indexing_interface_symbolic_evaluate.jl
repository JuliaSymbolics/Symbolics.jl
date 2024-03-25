using Symbolics
using SymbolicIndexingInterface
using Symbolics: Differential, Operator

@variables t x(t) y(t)
@variables p[1:3, 1:3] q[1:3]

bar(x, p) = p * x
@register_array_symbolic bar(x::AbstractVector, p::AbstractMatrix) begin
    size = size(x)
    eltype = promote_type(eltype(x), eltype(p))
end

D = Differential(t)

expr1 = x + y + D(x)
@test isequal(symbolic_evaluate(expr1, Dict(x => 3)), 3 + y + D(3))
@test isequal(symbolic_evaluate(expr1, Dict(x => 3); operator = Operator), 3 + y + D(x))
@test isequal(symbolic_evaluate(expr1, Dict(x => 1, D(x) => 2)), y + 3)
@test symbolic_evaluate(expr1, Dict(x => 1, D(x) => 2, y => 3)) == 6
@test isequal(symbolic_evaluate(expr1, Dict(x => 3, y => 3x), operator = Operator), 12 + D(x))
@test symbolic_evaluate(expr1, Dict(x => 3, y => 3x, D(x) => 2)) == 14

expr2 = bar(q, p)
@test isequal(symbolic_evaluate(expr2, Dict(p => ones(3, 3))), bar(q, ones(3, 3)))
@test symbolic_evaluate(expr2, Dict(p => ones(3, 3), q => ones(3))) == 3ones(3)

expr3 = bar(3q, 3p)
@test isequal(symbolic_evaluate(expr3, Dict(p => ones(3, 3))), bar(3q, 3ones(3, 3)))
@test symbolic_evaluate(expr3, Dict(p => ones(3, 3), q => ones(3))) == 27ones(3)

expr4 = D(x) ~ 3x + y
@test isequal(symbolic_evaluate(expr4, Dict(x => 3)), D(3) ~ 9 + y)
@test isequal(symbolic_evaluate(expr4, Dict(x => 3); operator = Operator), D(x) ~ y + 9)
@test isequal(symbolic_evaluate(expr4, Dict(x => 1, D(x) => 2)), 2 ~ 3 + y)
@test isequal(symbolic_evaluate(expr4, Dict(x => 1, D(x) => 2, y => 3)), 2 ~ 6)
