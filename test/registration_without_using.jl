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

using Symbolics, Test

# Both macros used to escape a bare `@wrapped` (and, before that, a bare
# `Symbolics.@wrapped`) into the caller. `wrap_func_expr` is now expanded inside
# Symbolics, so neither name has to be bound where the macro is invoked.
@testset "registration does not need `@wrapped` in the caller" begin
    for m in (RegistrationWithoutUsingImport, RegistrationWithoutUsingSelective)
        @test !isdefined(m, Symbol("@wrapped"))
    end
    @test !isdefined(RegistrationWithoutUsingSelective, :Symbolics)
end

# Registering into a bare `Module` is the strictest version of the above: only
# the macro itself is bound, so anything the expansion names by hand is missing.
@testset "registration into a module that only imports the macro" begin
    @variables x v[1:3]

    m = Module()
    Base.eval(m, :(using Symbolics: @register_symbolic))
    Base.eval(m, :(scalarfun(a, b, c) = a + b + c))
    # exercises the annotated-return-type path
    @test_nowarn Base.eval(m, :(@register_symbolic scalarfun(a::Real, b::Real, c::Real)::Real))
    @test Base.invokelatest(m.scalarfun, x, 1.0, 2.0) isa Num

    m2 = Module()
    Base.eval(m2, :(using Symbolics: @register_array_symbolic))
    Base.eval(m2, :(vand(x::AbstractVector) = x * x'))
    @test_nowarn Base.eval(m2, :(@register_array_symbolic vand(x::AbstractVector) begin
        size = (length(x), length(x))
        eltype = eltype(x)
    end))
    @test Base.invokelatest(m2.vand, v) isa Symbolics.Arr
end
