using Symbolics
import Symbolics: CallWithMetadata, getsource, getdefaultval, wrap, unwrap, getname
import SymbolicUtils: Term, symtype, FnType, BasicSymbolic, promote_symtype
using LinearAlgebra
using Test

@variables t
Symbolics.@register_symbolic fff(t)
@test isequal(fff(t), Symbolics.Num(Symbolics.Term{Real}(fff, [Symbolics.value(t)])))

const SymMatrix{T,N} =  Symmetric{T, AbstractArray{T, N}}
many_vars = @variables t=0 a=1 x[1:4]=2 y(t)[1:4]=3 w[1:4] = 1:4 z(t)[1:4] = 2:5 p(..)[1:4]

let
    @register_array_symbolic ggg(x::AbstractVector) begin
        container_type=SymMatrix
        size=(length(x) * 2, length(x) * 2)
        eltype=eltype(x)
    end false

    ## @variables

    gg = ggg(x)

    @test ndims(gg) == 2
    @test size(gg) == (8,8)
    @test eltype(gg) == Real
    @test symtype(unwrap(gg)) == SymMatrix{Real, 2}
    @test promote_symtype(ggg, symtype(unwrap(x))) == Any # no promote_symtype defined

    gg = ggg([a, 2a])
    @test ndims(gg) == 2
    @test size(gg) == (4, 4)
    @test eltype(gg) == Real
    @test symtype(unwrap(gg)) == SymMatrix{Real, 2}
    @test promote_symtype(ggg, Vector{symtype(typeof(a))}) == Any

    _a = unwrap(a)
    gg = ggg([_a, 2_a])
    @test ndims(gg) == 2
    @test size(gg) == (4, 4)
    @test eltype(gg) == Real
    @test symtype(unwrap(gg)) == SymMatrix{Real, 2}
    @test promote_symtype(ggg, Vector{symtype(typeof(a))}) == Any
end
let
    # redefine with promote_symtype
    @register_array_symbolic ggg(x::AbstractVector) begin
        container_type=SymMatrix
        size=(length(x) * 2, length(x) * 2)
        eltype=eltype(x)
    end
    @test promote_symtype(ggg, symtype(unwrap(x))) == SymMatrix{Real}
end

# ndims specified

# in terms of argument
@register_array_symbolic ggg(x::AbstractVector) begin
    container_type=SymMatrix
    size=(length(x) * 2, length(x) * 2)
    eltype=eltype(x)
    ndims = ndims(x) + 1
end
@test promote_symtype(ggg, symtype(unwrap(x))) == SymMatrix{Real, 2}

@register_array_symbolic ggg(x::AbstractVector) begin
    container_type=SymMatrix
    size=(length(x) * 2, length(x) * 2)
    eltype=eltype(x)
    ndims = 2
end
@test promote_symtype(ggg, symtype(unwrap(x))) == SymMatrix{Real, 2}

gg = ggg(x)

struct CanCallWithArray{T}
    params::T
end

ccwa = CanCallWithArray((length=10,))
@register_array_symbolic (c::CanCallWithArray)(x::AbstractArray, b::AbstractVector) begin
    size=(size(x, 1), length(b), c.params.length)
    eltype=Real
end false # without promote_symtype

hh = ccwa(gg, x)
@test size(hh) == (8,4,10)
@test eltype(hh) == Real
@test isequal(arguments(unwrap(hh)), unwrap.([gg, x]))

_args = [[a 2a; 4a 6a; 3a 5a], [4a, 6a]]
hh = ccwa(_args...)
@test size(hh) == (3, 2, 10)
@test eltype(hh) == Real
@test isequal(arguments(unwrap(hh)), unwrap.(_args))

@test all(t->getsource(t)[1] === :variables, many_vars)
@test getdefaultval(t) == 0
@test getdefaultval(a) == 1
@test getdefaultval(x) == 2
@test getdefaultval(x[1]) == 2
@test getdefaultval(y[2]) == 3
@test getdefaultval(w[2]) == 2
@test getdefaultval(w[4]) == 4
@test getdefaultval(z[3]) == 4

@test symtype(p) <: FnType{Tuple, Array{Real,1}}
@test promote_symtype(ccwa, symtype(unwrap(gg)), symtype(unwrap(x))) == Any
@test p(t)[1] isa Symbolics.Num


struct CanCallWithArray2{T}
    params::T
end

ccwa = CanCallWithArray2((length=10,))
@register_array_symbolic (c::CanCallWithArray2)(x::AbstractArray, b::AbstractVector) begin
    size=(size(x, 1), length(b), c.params.length)
    eltype=Real
end
@test promote_symtype(ccwa, symtype(unwrap(gg)), symtype(unwrap(x))) == AbstractArray{Real}

struct CanCallWithArray3{T}
    params::T
end

ccwa = CanCallWithArray3((length=10,))
# ndims specified
@register_array_symbolic (c::CanCallWithArray3)(x::AbstractArray, b::AbstractVector) begin
    size=(size(x, 1), length(b), c.params.length)
    eltype=Real
    ndims = 3
end
@test promote_symtype(ccwa, symtype(unwrap(gg)), symtype(unwrap(x))) == AbstractArray{Real, 3}

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
    vars = @variables t a b(a) c(..) x[1:2] y(t)[1:3] z(..)[1:2]
    vars2 = [Symbolics.rename(v, Symbol(Symbolics.getname(v), "_2")) for v in vars]

    for (v, v2) in zip(vars, vars2)
        @test typeof(v) == typeof(v2)
        @test Symbol(Symbolics.getname(v), "_2") == Symbolics.getname(v2)
    end


    @test Symbolics.getname(Symbolics.rename(y[2], :u)) === :u
end

let
    s = :y
    x = (1:2,1:3)
    t, y = @variables t $s(t)[x...]

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
    @test bar(unwrap(x), A()) isa BasicSymbolic
    @test typeof(baz(x, unwrap(y))) == Num
    @test typeof(baz(unwrap(x), unwrap(y))) <: BasicSymbolic
end

# 402#issuecomment-1074261734
Symbolics.@register_symbolic oof(x::AbstractVector)
Symbolics.@variables x[1:100]
@test oof(x) isa Num

# Register callable structs
# 806
struct B
    x::Float64
end
(b::B)(x) = b.x * x

Symbolics.@register_symbolic (b::B)(x)

@test B(2.0)(2.0) == 4.0
let
    foo = B(2.0)
    @variables x
    @test foo(x) isa Num
    @test foo(unwrap(x)) isa BasicSymbolic
end

@variables t y(t)
yy = Symbolics.variable(:y, T = Symbolics.FnType{Tuple{Any}, Real})
yyy = yy(t)
@test isequal(yyy, y)
@test yyy isa Num
@test y isa Num
yy = Symbolics.variable(:y, T = Symbolics.FnType{Tuple, Real})
yyy = yy(t)
@test !isequal(yyy, y)
@variables y(..)
@test isequal(yyy, y(t))

spam(x) = 2x
@register_symbolic spam(x::AbstractArray)

sym = spam([a, 2a])
@test sym isa Num
@test unwrap(sym) isa BasicSymbolic{Real}

fn_defaults = [print, min, max, identity, (+), (-), max, sum, vcat, (*)]
fn_names = [Symbol(:f, i) for i in 1:10]

struct VariableFoo end
Symbolics.option_to_metadata_type(::Val{:foo}) = VariableFoo

function test_all_functions(fns)
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10 = fns
    @variables x y::Int z::Function w[1:3, 1:3] v[1:3, 1:3]::String
    @test f1 isa CallWithMetadata{FnType{Tuple, Real}}
    @test all(x -> symtype(x) <: Real, [f1(), f1(1), f1(x), f1(x, y), f1(x, y, x+y)])
    @test f2 isa CallWithMetadata{FnType{Tuple{Any, Vararg}, Int}}
    @test all(x -> symtype(x) <: Int, [f2(1), f2(z), f2(x), f2(x, y), f2(x, y, x+y)])
    @test_throws ErrorException f2()
    @test f3 isa CallWithMetadata{FnType{Tuple, Real, typeof(max)}}
    @test all(x -> symtype(x) <: Real, [f3(), f3(1), f3(x), f3(x, y), f3(x, y, x+y)])
    @test f4 isa CallWithMetadata{FnType{Tuple{Int}, Real}}
    @test all(x -> symtype(x) <: Real, [f4(1), f4(y), f4(2y)])
    @test_throws ErrorException f4(x)
    @test f5 isa CallWithMetadata{FnType{Tuple{Int, Vararg{Int}}, Real}}
    @test all(x -> symtype(x) <: Real, [f5(1), f5(y), f5(y, y), f5(2, 3)])
    @test_throws ErrorException f5(x)
    @test f6 isa CallWithMetadata{FnType{Tuple{Int, Int}, Int}}
    @test all(x -> symtype(x) <: Int, [f6(1, 1), f6(y, y), f6(1, y), f6(y, 1)])
    @test_throws ErrorException f6()
    @test_throws ErrorException f6(1)
    @test_throws ErrorException f6(x, y)
    @test_throws ErrorException f6(y)
    @test f7 isa CallWithMetadata{FnType{Tuple{Int, Int}, Int, typeof(max)}}
    # call behavior tested by f6
    @test f8 isa CallWithMetadata{FnType{Tuple{Function, Vararg}, Real, typeof(sum)}}
    @test all(x -> symtype(x) <: Real, [f8(z), f8(z, x), f8(identity), f8(identity, x)])
    @test_throws ErrorException f8(x)
    @test_throws ErrorException f8(1)
    @test f9 isa CallWithMetadata{FnType{Tuple, Vector{Real}}}
    @test all(x -> symtype(unwrap(x)) <: Vector{Real} && size(x) == (3,), [f9(), f9(1), f9(x), f9(x + y), f9(z), f9(1, x)])
    @test f10 isa CallWithMetadata{FnType{Tuple{Matrix{<:Real}, Matrix{<:Real}}, Matrix{Real}, typeof(*)}}
    @test all(x -> symtype(unwrap(x)) <: Matrix{Real} && size(x) == (3, 3), [f10(w, w), f10(w, ones(3, 3)), f10(ones(3, 3), ones(3, 3)), f10(w + w, w)])
    @test_throws ErrorException f10(w, v)
end

function test_functions_defaults(fns)
    for (fn, def) in zip(fns, fn_defaults)
        @test Symbolics.getdefaultval(fn, nothing) == def
    end
end

function test_functions_metadata(fns)
    for (i, fn) in enumerate(fns)
        @test Symbolics.getmetadata(fn, VariableFoo, nothing) == i
    end
end

fns = @test_nowarn @variables begin
    f1(..)
    f2(::Any, ..)::Int
    (f3::typeof(max))(..)
    f4(::Int)
    f5(::Int, (..)::Int)
    f6(::Int, ::Int)::Int
    (f7::typeof(max))(::Int, ::Int)::Int
    (f8::typeof(sum))(::Function, ..)
    f9(..)[1:3]
    (f10::typeof(*))(::Matrix{<:Real}, ::Matrix{<:Real})[1:3, 1:3]
    # f11[1:3](::Int)::Int
end

test_all_functions(fns)

fns = @test_nowarn @variables begin
    f1(..) = fn_defaults[1]
    f2(::Any, ..)::Int = fn_defaults[2]
    (f3::typeof(max))(..) = fn_defaults[3]
    f4(::Int) = fn_defaults[4]
    f5(::Int, (..)::Int) = fn_defaults[5]
    f6(::Int, ::Int)::Int = fn_defaults[6]
    (f7::typeof(max))(::Int, ::Int)::Int = fn_defaults[7]
    (f8::typeof(sum))(::Function, ..) = fn_defaults[8]
    f9(..)[1:3] = fn_defaults[9]
    (f10::typeof(*))(::Matrix{<:Real}, ::Matrix{<:Real})[1:3, 1:3] = fn_defaults[10]
end

test_all_functions(fns)
test_functions_defaults(fns)

fns = @variables begin
    f1(..) = fn_defaults[1], [foo = 1]
    f2(::Any, ..)::Int = fn_defaults[2], [foo = 2;]
    (f3::typeof(max))(..) = fn_defaults[3], [foo = 3;]
    f4(::Int) = fn_defaults[4], [foo = 4;]
    f5(::Int, (..)::Int) = fn_defaults[5], [foo = 5;]
    f6(::Int, ::Int)::Int = fn_defaults[6], [foo = 6;]
    (f7::typeof(max))(::Int, ::Int)::Int = fn_defaults[7], [foo = 7;]
    (f8::typeof(sum))(::Function, ..) = fn_defaults[8], [foo = 8;]
    f9(..)[1:3] = fn_defaults[9], [foo = 9;]
    (f10::typeof(*))(::Matrix{<:Real}, ::Matrix{<:Real})[1:3, 1:3] = fn_defaults[10], [foo = 10;]
end

test_all_functions(fns)
test_functions_defaults(fns)
test_functions_metadata(fns)

fns = @test_nowarn @variables begin
    f1(..), [foo = 1,]
    f2(::Any, ..)::Int, [foo = 2,]
    (f3::typeof(max))(..), [foo = 3,]
    f4(::Int), [foo = 4,]
    f5(::Int, (..)::Int), [foo = 5,]
    f6(::Int, ::Int)::Int, [foo = 6,]
    (f7::typeof(max))(::Int, ::Int)::Int, [foo = 7,]
    (f8::typeof(sum))(::Function, ..), [foo = 8,]
    f9(..)[1:3], [foo = 9,]
    (f10::typeof(*))(::Matrix{<:Real}, ::Matrix{<:Real})[1:3, 1:3], [foo = 10,]
end

test_all_functions(fns)
test_functions_metadata(fns)

fns = @test_nowarn @variables begin
    $(fn_names[1])(..)
    $(fn_names[2])(::Any, ..)::Int
    ($(fn_names[3])::typeof(max))(..)
    $(fn_names[4])(::Int)
    $(fn_names[5])(::Int, (..)::Int)
    $(fn_names[6])(::Int, ::Int)::Int
    ($(fn_names[7])::typeof(max))(::Int, ::Int)::Int
    ($(fn_names[8])::typeof(sum))(::Function, ..)
    $(fn_names[9])(..)[1:3]
    ($(fn_names[10])::typeof(*))(::Matrix{<:Real}, ::Matrix{<:Real})[1:3, 1:3]
end

test_all_functions(fns)
