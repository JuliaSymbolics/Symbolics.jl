module RegistrationWithoutUsingImport
import Symbolics
Symbolics.@register_symbolic foo(x, y)
Symbolics.@register_array_symbolic bar(x::AbstractVector) begin
    size = (length(x), length(x))
end false
end

module RegistrationWithoutUsingSelective
import Symbolics: @register_array_symbolic, @register_symbolic
@register_symbolic goo(x, y)
@register_array_symbolic hoo(x::AbstractVector) begin
    size = (length(x), length(x))
end false
end
