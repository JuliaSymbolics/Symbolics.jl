import Symbolics
Symbolics.@register_symbolic foo(x, y)
Symbolics.@register_array_symbolic bar(x::AbstractVector) begin
    size = (length(x), length(x))
end false
