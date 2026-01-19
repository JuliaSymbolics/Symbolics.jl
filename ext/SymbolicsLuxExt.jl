module SymbolicsLuxExt

using Lux
using Symbolics
using Lux.LuxCore
using Lux.Random: AbstractRNG, default_rng
using Symbolics: SymbolicT, VartypeT
using SymbolicUtils
using SymbolicUtils: symtype, unwrap

@static if isdefined(Lux.NilSizePropagation, :recursively_nillify)
    function Lux.NilSizePropagation.recursively_nillify(x::SymbolicUtils.BasicSymbolic)
        @assert SymbolicUtils.symtype(x) <: Vector{<:Real}
        Lux.NilSizePropagation.recursively_nillify(Symbolics.Arr{Num, 1}(x))
    end
end

# SymbolicUtils@4 requires implementing `promote_symtype(::typeof(f), sh1::Shapet, sh2::ShapeT, ...)`
# for `substitute` and other functions to work correctly. Effectively, it requires being
# able to compute the shape of the output from the function and shapes of the inputs. With
# `LuxCore.stateless_apply`, this becomes difficult since it requires using `LuxCore.outputsize`,
# which needs the model as input. The code below takes a similar approach to the
# implementation of `map` and `mapreduce` in SymbolicUtils. A functor
# `LuxStatelessApplicator` is used as the operation, with `x` and `ps` as the arguments.
# The functor contains the (possibly symbolic) model and verifies that a concrete model
# can be obtained from it. This allows implementing `promote_shape`, since the model
# is available as part of the function.
#
# `(::SymbolicUtils.Substituter)` is implemented to be able to substitute the model
# inside `LuxStatelessApplicator`. `SymbolicUtils.Code.function_to_expr` is implemented
# so that codegen generates a call to `LuxCore.stateless_apply`. Calling the functor
# itself also calls `LuxCore.stateless_apply`, which is necessary for `substitute` with
# `fold = Val(true)` to work correctly. In case some arguments to the functor are
# symbolic, it will again hit the manual registration dispatch and create a symbolic
# expression representing the model forward pass.

struct LuxStatelessApplicator
    model::Union{LuxCore.AbstractLuxLayer, SymbolicT}

    LuxStatelessApplicator(model::LuxCore.AbstractLuxLayer) = new(model)
    function LuxStatelessApplicator(model::SymbolicT)
        if SymbolicUtils.isconst(model)
            return LuxStatelessApplicator(unwrap_const(model)::LuxCore.AbstractLuxLayer)
        end
        @assert symtype(model) <: LuxCore.AbstractLuxLayer
        @assert Symbolics.getdefaultval(model, nothing) isa LuxCore.AbstractLuxLayer
        new(model)
    end
end

Base.nameof(::LuxStatelessApplicator) = :LuxStatelessApplicator

function __get_model(f::LuxStatelessApplicator)
    sym = f.model
    if sym isa LuxCore.AbstractLuxLayer
        return sym
    else
        return Symbolics.getdefaultval(sym::SymbolicT)::LuxCore.AbstractLuxLayer
    end
end

function SymbolicUtils.Code.function_to_expr(op::LuxStatelessApplicator, O, st)
    rw = get(st.rewrites, O, nothing)
    rw === nothing || return rw

    return Expr(:call, LuxCore.stateless_apply, SymbolicUtils.Code.toexpr(op.model, st),
                SymbolicUtils.Code.toexpr(O.args[1], st), SymbolicUtils.Code.toexpr(O.args[2], st))
end

function (s::SymbolicUtils.Substituter)(f::LuxStatelessApplicator)
    return LuxStatelessApplicator(s(f.model))
end

function (f::LuxStatelessApplicator)(x, ps)
    return LuxCore.stateless_apply(f.model, x, ps)
end

for modelT in [LuxCore.AbstractLuxLayer, SymbolicT]
    @eval function LuxCore.stateless_apply(model::$modelT, x::Union{SymbolicT, AbstractArray{SymbolicT}, AbstractArray{Num}}, ps)
        f = LuxStatelessApplicator(model)
        x = SymbolicUtils.BSImpl.Const{VartypeT}(x)
        ps = SymbolicUtils.BSImpl.Const{VartypeT}(ps)
        T = SymbolicUtils.promote_symtype(f, symtype(x), symtype(ps))
        sh = SymbolicUtils.promote_shape(f, SymbolicUtils.shape(x), SymbolicUtils.shape(ps))
        args = SymbolicUtils.ArgsT{VartypeT}((x, ps))
        return SymbolicUtils.BSImpl.Term{VartypeT}(f, args; type = T, shape = sh)
    end

    @eval function LuxCore.stateless_apply(model::$modelT, x::Symbolics.Arr{T, N}, ps) where {T, N}
        Symbolics.Arr{T, N}(LuxCore.stateless_apply(model, unwrap(x), ps))
    end
end

function SymbolicUtils.promote_symtype(::LuxStatelessApplicator, x::SymbolicUtils.TypeT, ps::SymbolicUtils.TypeT)
    @assert x <: AbstractArray
    if x <: Array
        return Array{Real, x.parameters[2]::Int}
    else
        return Array{Real, ndims(x)::Int}
    end
end

function SymbolicUtils.promote_shape(f::LuxStatelessApplicator, x::SymbolicUtils.ShapeT, ::SymbolicUtils.ShapeT)
    @nospecialize x ps
    model = __get_model(f)
    x = x::SymbolicUtils.ShapeVecT
    mockup = SymbolicUtils.Sym{Symbolics.VartypeT}(:__Lux_mockup_x; shape = x, type = Array{Real, length(x)})
    mockup = Symbolics.Arr{Num, length(x)}(mockup)
    sz = LuxCore.outputsize(model, mockup, default_rng())
    sh = SymbolicUtils.ShapeVecT()
    for i in sz
        push!(sh, 1:i)
    end
    return sh
end

end
