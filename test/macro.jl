using Symbolics
import Symbolics: getsource, getdefaultval, wrap, unwrap, getname
import SymbolicUtils: Term, symtype, FnType
using Test

@variables t
Symbolics.@register_symbolic fff(t)
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


let
    vars = @variables t a b(a) c(..) x[1:2] y[1:3](t) z[1:2](..)
    vars2 = [Symbolics.rename(v, Symbol(Symbolics.getname(v), "_2")) for v in vars]

    for (v, v2) in zip(vars, vars2)
        @test typeof(v) == typeof(v2)
        @test Symbol(Symbolics.getname(v), "_2") == Symbolics.getname(v2)
    end

    @test Symbolics.getname(getmetadata(z[2](t), Symbolics.GetindexParent)) === :z
    @test Symbolics.getname(getmetadata(vars2[end][2](t), Symbolics.GetindexParent)) === :z_2

    @test Symbolics.getname(Symbolics.rename(y[2], :u)) === :u
end

let
    s = :y
    x = (1:2,1:3)
    t, y = @variables t $s[x...](t)

    @test ndims(y) == 2
    @test size(y) == (2,3)
end


# Edge case when no symbolic values are passed
struct A end
Symbolics.@register_symbolic bar(t, x::A)
Symbolics.@register_symbolic baz(x, y)
if !@isdefined(bar_catchall_defined)
    @test_throws MethodError bar(0.1, A())
    @test_throws MethodError bar(Num(0.1), A())
else
    @warn("skipping 2 tests because this file was run more than once")
end

bar_catchall_defined = true
bar(t, x::A) = 1
@test bar(1, A()) == 1

let
    @variables x y
    @test bar(x, A()) isa Num
    @test bar(unwrap(x), A()) isa Term
    @test typeof(baz(x, unwrap(y))) == Num
    @test typeof(baz(unwrap(x), unwrap(y))) <: Term
end

# 402#issuecomment-1074261734
Symbolics.@register_symbolic oof(x::AbstractVector)
Symbolics.@variables x[1:100]
@test oof(x) isa Num
