module SymbolicsLatexifyExt

using Symbolics
using Latexify
using LaTeXStrings
using TermInterface
using SymbolicUtils
using Symbolics: value, hide_lhs, wrap
using MacroTools: postwalk, prewalk
using SymbolicUtils: BSImpl, FnType, unwrap, symtype, BasicSymbolic
using Moshi.Match: @match

# metadata to specify how to format syms
struct SymLatexWrapper end
Symbolics.option_to_metadata_type(::Val{:latexwrapper}) = SymLatexWrapper

function default_latex_wrapper(sym)
    length(sym) <= 1 && return sym
    return string("\\mathtt{", sym, "}")
end

prettify_expr(expr) = expr
prettify_expr(f::Function) = nameof(f)
prettify_expr(expr::Expr) = Expr(expr.head, prettify_expr.(expr.args)...)

function cleanup_exprs(ex)
    return postwalk(x -> iscall(x) && length(arguments(x)) == 0 ? operation(x) : x, ex)
end

# Keep custom-wrapper call arguments as Expr nodes so the outer Latexify
# traversal applies recipe/caller options (index, fmt, mult_symbol, snakecase).
# Always `:block`-wrap Expr args: `:latexifymerge` wraps any non-`:none` child
# in an extra `\left(...\right)`, and we already supply the call parentheses.
function _latexstring_call_to_merge(ex::Expr)
    name = ex.args[1]
    length(ex.args) == 1 && return name
    body = Expr(:latexifymerge, name, "\\left( ")
    for (i, a) in enumerate(ex.args[2:end])
        i > 1 && (body = Expr(:latexifymerge, body, ", "))
        child = a isa Expr ? Expr(:block, a) : a
        body = Expr(:latexifymerge, body, child)
    end
    return Expr(:latexifymerge, body, " \\right)")
end

function latexify_derivatives(ex)
    # Latexify does not parenthesize `^` when the base is a LaTeXString-headed
    # call (`_getoperation` only recognizes Symbol heads). Mark those first.
    ex = prewalk(ex) do x
        if Meta.isexpr(x, :call) && x.args[1] == :^ && length(x.args) >= 3
            base = x.args[2]
            if Meta.isexpr(base, :call) && base.args[1] isa LaTeXString && length(base.args) > 1
                return Expr(:call, :^, Expr(:call, :_latexfenced, base), x.args[3])
            end
        end
        return x
    end
    # Pass 1: derivatives/integrals while unary custom calls are still `:call`
    # nodes, so `D(f(x))` keeps `f(x)` in the fraction numerator (needed when
    # the derivative is later multiplied by another factor).
    ex = postwalk(ex) do x
        Meta.isexpr(x, :call) || return x
        if x.args[1] == :_derivative
            num, den, deg = x.args[2:end]
            dsym = "\\mathrm{d}$(deg == 1 ? "" : "^{$deg}")"
            den_ls = diffdenom(den)
            if Meta.isexpr(num, :call) && length(num.args) == 2 && num.args[1] !== :*
                return Expr(:call, :/, Expr(:latexifymerge, dsym, _latexify_merge_child(num)), den_ls)
            else
                return Expr(
                    :latexifymerge,
                    LaTeXString("\\frac{$dsym}{$(den_ls.s)} ~ "),
                    _latexify_merge_child(num)
                )
            end
        elseif x.args[1] === :_integral
            lower, upper, var_of_int, integrand = x.args[2:end]
            body = Expr(:latexifymerge, "\\int_{", _latexify_merge_child(lower))
            body = Expr(:latexifymerge, body, Expr(:latexifymerge, "}^{", _latexify_merge_child(upper)))
            body = Expr(:latexifymerge, body, "} ~ ")
            body = Expr(:latexifymerge, body, _latexify_merge_child(var_of_int))
            body = Expr(:latexifymerge, body, Expr(:latexifymerge, " ~ ", _latexify_merge_child(integrand)))
            return body
        elseif x.args[1] == :^ && length(x.args) == 3 && _latexify_power_base_needs_parens(x.args[2])
            # `:latexifymerge` has no precedence; parenthesise a differential/integral
            # form used as a power base.
            return Expr(
                :call, :^,
                Expr(:latexifymerge, "\\left( ", x.args[2], " \\right)"),
                x.args[3]
            )
        elseif x.args[1] === :_textbf
            ls = latexify(latexify_derivatives(sorted_arguments(x)[1])).s
            return "\\textbf{" * strip(ls, '\$') * "}"
        else
            return x
        end
    end
    # Pass 2: convert custom-wrapper calls to `:latexifymerge` and apply fences.
    return postwalk(ex) do x
        Meta.isexpr(x, :call) || return x
        if x.args[1] === :_latexfenced
            inner = x.args[2]
            if Meta.isexpr(inner, :call) && inner.args[1] isa LaTeXString
                inner = _latexstring_call_to_merge(inner)
            end
            return Expr(:latexifymerge, "\\left( ", inner, " \\right)")
        elseif x.args[1] isa LaTeXString
            return _latexstring_call_to_merge(x)
        else
            return x
        end
    end
end

# `:latexifymerge` parenthesises any child with a non-`:none` operation. Leave
# binary `+`/`*`/`/`/`-` bare so they stay grouped; wrap other `Expr` children
# in `:block` so calls, refs and powers stay bare. Non-`Expr` atoms are already
# `:none` and must not be block-wrapped (Latexify treats the block arg as `op`).
function _latexify_needs_merge_parens(ex)
    Meta.isexpr(ex, :call) || return false
    op = ex.args[1]
    (op isa Symbol && Base.isoperator(op)) || return false
    op === :^ && return false
    return length(ex.args) >= 3
end

function _latexify_merge_child(ex)
    _latexify_needs_merge_parens(ex) && return ex
    return ex isa Expr ? Expr(:block, ex) : ex
end

function _latexify_power_base_needs_parens(base)
    Meta.isexpr(base, :latexifymerge) || return false
    a1 = base.args[1]
    if a1 isa AbstractString || a1 isa LaTeXString
        s = a1 isa LaTeXString ? a1.s : a1
        return startswith(s, "\\frac") || startswith(s, "\\int")
    end
    return _latexify_power_base_needs_parens(a1)
end

# `latexify_derivatives` can collapse a top-level node into a bare `String` (e.g.
# `_textbf(...)` -> "\\textbf{...}", reachable both from the array-symbol branch and from
# a custom `_toexpr_metadata`/`_toexpr_op` hook). Latexify would try to re-parse such a
# `String` as an expression and fail, so wrap a top-level string as a `LaTeXString`, which
# is emitted verbatim. Strings nested inside an `Expr` are left untouched.
_as_latexstring(x) = x
_as_latexstring(x::AbstractString) = LaTeXString(x)

recipe(n) = _as_latexstring(latexify_derivatives(cleanup_exprs(_toexpr(n))))

function align_side(n)
    wrapped = wrap(n)
    return wrapped isa Symbolics.Arr ? recipe(n) : wrapped
end

@latexrecipe function f(n::Num)
    env --> :equation
    mult_symbol --> "~"
    fmt --> FancyNumberFormatter(5)
    index --> :subscript
    snakecase --> true
    safescripts --> true

    return recipe(value(n))
end

@latexrecipe function f(z::Complex{Num})
    env --> :equation
    mult_symbol --> "~"
    index --> :subscript

    iszero(z.im) && return :($(recipe(value(z.re))))
    iszero(z.re) && return :($(recipe(value(z.im))) * $im)
    return :($(recipe(value(z.re))) + $(recipe(value(z.im))) * $im)
end

@latexrecipe function f(n::Function)
    env --> :equation
    mult_symbol --> "~"
    index --> :subscript

    return nameof(n)
end


@latexrecipe function f(n::Symbolics.Arr)
    env --> :equation
    mult_symbol --> "~"
    index --> :subscript

    return value(n)
end

@latexrecipe function f(n::Symbolics.CallAndWrap)
    env --> :equation
    mult_symbol --> "~"
    index --> :subscript

    return n.f
end

@latexrecipe function f(n::SymbolicUtils.BasicSymbolic)
    env --> :equation
    mult_symbol --> "~"
    index --> :subscript

    return recipe(n)
end

@latexrecipe function f(eqs::Vector{Equation})
    index --> :subscript
    has_connections = any(x -> hide_lhs(value(x.lhs)), eqs)
    if has_connections
        env --> :equation
        return map(first ∘ first ∘ Latexify.apply_recipe, eqs)
    else
        env --> :align
        return align_side.(getfield.(eqs, :lhs)), align_side.(getfield.(eqs, :rhs))
    end
end

@latexrecipe function f(eq::Equation)
    env --> :equation
    index --> :subscript

    if hide_lhs(value(eq.lhs)) || !(value(eq.lhs) isa Union{Number, AbstractArray, BasicSymbolic})
        return value(eq.rhs)
    else
        return Expr(:(=), recipe(eq.lhs), recipe(eq.rhs))
    end
end

Base.show(io::IO, ::MIME"text/latex", x::Symbolics.RCNum) = print(io, "\$\$ " * latexify(x) * " \$\$")
Base.show(io::IO, ::MIME"text/latex", x::SymbolicUtils.BasicSymbolic) = print(io, "\$\$ " * latexify(x) * " \$\$")
Base.show(io::IO, ::MIME"text/latex", x::Equation) = print(io, "\$\$ " * latexify(x) * " \$\$")
Base.show(io::IO, ::MIME"text/latex", x::Vector{Equation}) = print(io, "\$\$ " * latexify(x) * " \$\$")
Base.show(io::IO, ::MIME"text/latex", x::AbstractArray{<:Symbolics.RCNum}) = print(io, "\$\$ " * latexify(x) * " \$\$")

# Iterate a node's metadata, dispatching to the `Symbolics._toexpr_metadata` hook for
# each context. The per-context hook (and the `Symbolics._toexpr_op` hook used below) is
# defined in `Symbolics` so downstream packages can extend it via `import Symbolics`
# without reaching into this extension with `Base.get_extension`.
function Symbolics._toexpr_metadata(O; latexwrapper = default_latex_wrapper)
    md = SymbolicUtils.metadata(O)
    md isa AbstractDict || return nothing
    for (ctx, val) in md
        out = Symbolics._toexpr_metadata(O, ctx, val; latexwrapper)
        out === nothing || return out
    end
    return nothing
end

# `_toexpr` is only used for latexify
function _toexpr(O; latexwrapper = default_latex_wrapper)
    O = unwrap(O)
    SymbolicUtils.isconst(O) && return value(O)
    custom = Symbolics._toexpr_metadata(O; latexwrapper)
    custom === nothing || return custom
    return _toexpr_plain(O; latexwrapper)
end

function _toexpr_plain(O; latexwrapper = default_latex_wrapper)
    if SymbolicUtils.ismul(O)
        m = O
        numer = Any[]
        denom = Any[]

        # We need to iterate over each term in m, ignoring the numeric coefficient.
        # This iteration needs to be stable, so we can't iterate over m.dict.
        for term in Iterators.drop(sorted_arguments(m), isone(m.coeff) ? 0 : 1)
            if !SymbolicUtils.ispow(term)
                push!(numer, _toexpr(term))
                continue
            end
            base, pow = arguments(term)
            pow = value(pow)
            isneg = (pow isa Number && pow < 0) || (iscall(pow) && operation(pow) === (-) && length(arguments(pow)) == 1)
            if !isneg
                if SymbolicUtils._isone(pow)
                    pushfirst!(numer, _toexpr(base))
                else
                    pushfirst!(numer, Expr(:call, :^, _toexpr(base), _toexpr(pow)))
                end
            else
                newpow = -1 * pow
                if SymbolicUtils._isone(newpow)
                    pushfirst!(denom, _toexpr(base))
                else
                    pushfirst!(denom, Expr(:call, :^, _toexpr(base), _toexpr(newpow)))
                end
            end
        end

        if !isreal(m.coeff)
            numer_expr = Expr(:call, :*, m.coeff, numer...)
        elseif isempty(numer) || !isone(abs(m.coeff))
            numer_expr = Expr(:call, :*, abs(m.coeff), numer...)
        else
            numer_expr = length(numer) > 1 ? Expr(:call, :*, numer...) : numer[1]
        end

        if isempty(denom)
            frac_expr = numer_expr
        else
            denom_expr = length(denom) > 1 ? Expr(:call, :*, denom...) : denom[1]
            frac_expr = Expr(:call, :/, numer_expr, denom_expr)
        end

        if isreal(m.coeff) && real(m.coeff) < 0
            return Expr(:call, :-, frac_expr)
        else
            return frac_expr
        end
    end
    if SymbolicUtils.issym(O)
        sym = string(nameof(O))
        sym = replace(sym, Symbolics.NAMESPACE_SEPARATOR => ".")

        # override if the sym has its own latex wrapper
        has_custom_wrapper = hasmetadata(O, SymLatexWrapper)
        symwrapper = has_custom_wrapper ? getmetadata(O, SymLatexWrapper) : latexwrapper
        sym = symwrapper(sym)
        # Custom wrappers supply raw LaTeX; emit LaTeXString so snakecase=true does not
        # escape `_`/`^`. Keep Symbol for the default wrapper so unannotated names match.
        if has_custom_wrapper || latexwrapper !== default_latex_wrapper
            return LaTeXString(sym)
        else
            return Symbol(sym)
        end
    end
    !iscall(O) && return O

    op = operation(O)
    args = sorted_arguments(O)
    latexwrapper = hasmetadata(O, SymLatexWrapper) ? getmetadata(O, SymLatexWrapper) :
        default_latex_wrapper

    custom_op = Symbolics._toexpr_op(op, args; latexwrapper)
    custom_op === nothing || return custom_op

    if (op === (*)) && (args[1] === -1)
        arg_mul = Expr(:call, :(*), _toexpr(args[2:end])...)
        return Expr(:call, :(-), arg_mul)
    end

    if op isa Differential
        #  DERIVATIVES LOGIC
        num = args[1]
        diff_var = op.x

        deg = op.order

        while iscall(num) && operation(num) isa Differential && isequal(operation(num).x, diff_var)
            inner_op = operation(num)
            deg += inner_op.order
            num = arguments(num)[1]
        end

        den = deg > 1 ? (diff_var^deg) : diff_var
        return :(_derivative($(_toexpr(num)), $den, $deg))

    elseif op isa Integral
        lower = op.domain.domain.left
        upper = op.domain.domain.right
        vars = op.domain.variables
        integrand = args[1]
        var = if vars isa Tuple
            Expr(:call, :(*), _toexpr(vars...))
        else
            _toexpr(vars)
        end
        return Expr(:call, :_integral, _toexpr(lower), _toexpr(upper), vars, _toexpr(integrand))
    elseif symtype(op) <: FnType
        isempty(args) && return nameof(op)
        return Expr(:call, _toexpr(op; latexwrapper), _toexpr(args)...)
    elseif op === getindex && symtype(args[1]) <: AbstractArray
        return getindex_to_symbol(O)
    elseif op === (\)
        return :(solve($(_toexpr(args[1])), $(_toexpr(args[2]))))
    elseif SymbolicUtils.issym(op) && SymbolicUtils.symtype(op) <: AbstractArray
        return :(_textbf($(nameof(op))))
    elseif op === identity
        return _toexpr(only(args)) # suppress identity transformations (e.g. "identity(π)" -> "π")
    end
    return Expr(:call, Symbol(op), _toexpr(args; latexwrapper)...)
end
_toexpr(x::Integer; latexwrapper = default_latex_wrapper) = x
_toexpr(x::AbstractFloat; latexwrapper = default_latex_wrapper) = x

function _toexpr(eq::Equation; latexwrapper = default_latex_wrapper)
    return Expr(:(=), _toexpr(eq.lhs), _toexpr(eq.rhs))
end

_toexpr(eqs::AbstractArray; latexwrapper = default_latex_wrapper) = map(eq -> _toexpr(eq), eqs)
_toexpr(x::Num; latexwrapper = default_latex_wrapper) = _toexpr(value(x))

function getindex_to_symbol(t)
    @assert iscall(t) && operation(t) === getindex && SymbolicUtils.symtype(sorted_arguments(t)[1]) <: AbstractArray
    args = sorted_arguments(t)
    idxs = args[2:end]
    O = args[1]
    latexwrapper = (O isa SymbolicUtils.BasicSymbolic && hasmetadata(O, SymLatexWrapper)) ? getmetadata(O, SymLatexWrapper) :
        default_latex_wrapper

    # this is to ensure X(t)[1] becomes X_1(t) in Latex
    if iscall(O) && SymbolicUtils.issym(operation(O))
        oop = operation(O)
        oargs = sorted_arguments(O)
        return :($(_toexpr(oop; latexwrapper))[$(idxs...)]($(_toexpr(oargs)...)))
    else
        return :($(_toexpr(O; latexwrapper))[$(idxs...)])
    end
end

function diffdenom(e)
    e = unwrap(e)
    return if SymbolicUtils.issym(e)
        LaTeXString("\\mathrm{d}$e")
    elseif SymbolicUtils.ispow(e)
        base, expo = arguments(e)
        suffix = SymbolicUtils._isone(expo) ? "" : "^{$(expo)}"
        LaTeXString("\\mathrm{d}$(base)$(suffix)")
    elseif SymbolicUtils.ismul(e)
        LaTeXString(prod(diffdenom(arg).s for arg in arguments(e)))
    else
        LaTeXString("\\mathrm{d}$e")
    end
end

end
