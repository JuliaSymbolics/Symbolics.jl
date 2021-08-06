using Symbolics
import Symbolics: getsource, getdefaultval, wrap, unwrap, getname
import SymbolicUtils: Term, symtype, FnType
using Test

@variables t
Symbolics.@register fff(t)
@test isequal(fff(t), Symbolics.Num(Symbolics.Term{Real}(fff, [Symbolics.value(t)])))

## @variables

many_vars = @variables t=0 a=1 x[1:4]=2 y[1:4](t)=3 w[1:4] = 1:4 z[1:4](t) = 2:5 p[1:4](..)

@test all(t->getsource(t)[1] === :variables, many_vars)
@test getdefaultval(t) == 0
@test getdefaultval(a) == 1
@test_throws ErrorException getdefaultval(x)
@test getdefaultval(x[1]) == 2
@test getdefaultval(y[2]) == 3
@test getdefaultval(w[2]) == 2
@test getdefaultval(w[4]) == 4
@test getdefaultval(z[3]) == 4

nxt = Namespace(unwrap(x), unwrap(t))
nxt_1 = Namespace(nxt, unwrap(x))
@test getname(nxt) == Symbol(:x, :(.), :t)
@test getname(nxt_1) == Symbol(:x, :(.), :t, :(.), :x)

@test p[1] isa Symbolics.CallWithMetadata
@test symtype(p[1]) <: FnType{Tuple, Real}
@test p[1](t) isa Symbolics.Num


## Wrapper types

abstract type AbstractFoo{T} end

struct Foo{T}<:AbstractFoo{T}
end
@symbolic_wrap struct FooWrap{T} <: AbstractFoo{T}
    val::Foo{T}
end

Symbolics.unwrap(r::FooWrap) = r.val

@test Symbolics.wrapper_type(AbstractFoo) == FooWrap
@test Symbolics.wrapper_type(AbstractFoo{Int}) == FooWrap # removes the type param.
                                                          # If required,
                                                          # a new method must be defined
                                                          # by the user of @symbolic_wrap

x = Foo{Int}()
@test wrap(x) isa FooWrap{Int}
@test unwrap(wrap(x)) isa Foo{Int}

@wrapped function foo(f::AbstractFoo, x::Real)
    x > 1 ? (x > 5 ? "hi" : Foo{Int}()) : 0
end

@test !applicable(foo, x, 2)
@test applicable(foo, x, wrap(2))
@test applicable(foo, wrap(x), 2)
@test applicable(foo, wrap(x), wrap(2))

@test foo(x, wrap(2)) isa FooWrap
@test foo(x, wrap(1)) isa Num
@test foo(x, wrap(6)) isa String
