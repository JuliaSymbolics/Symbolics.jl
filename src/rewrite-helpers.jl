"""
    replace(expr::Symbolic, rules...)
Walk the expression and replace subexpressions according to `rules`. `rules`
could be rules constructed with `@rule`, a function, or a pair where the
left hand side is matched with equality (using `isequal`) and is replaced by the right hand side.

Rules will be applied left-to-right simultaneously,
so only one pattern will be applied to any subexpression,
and the patterns will only be applied to the input text,
not the replacements.

Set `fixpoint = true` to repeatedly apply rules until no
change to the expression remains to be made.
"""
function _replace(expr::Symbolic, rules...; fixpoint=false)
    rs = map(r -> r isa Pair ? (x -> isequal(x, r[1]) ? r[2] : nothing) : r, rules)
    R = Prewalk(Chain(rs))
    if fixpoint
        Fixpoint(R)(expr)
    else
        R(expr)
    end
end
# Fix ambiguity
function Base.replace(expr::Num, r::Pair, rules::Pair...)
    _replace(unwrap(expr), r, rules...)
end

function Base.replace(expr::Num, rules...)
    _replace(unwrap(expr), rules...)
end

function Base.replace(expr::Symbolic, r, rules...)
    _replace(expr, r, rules)
end

"""
    occursin(x, y)
Checks whether `x` occurs in `y`. Parses `y`, returning `true` on any occurrence of x.

Example:
```julia
@syms x y z
occursin(x, y + z*(3+x)) # returns `true`.
occursin(x, log(y) + 10*(z-y) # returns `false`.
```
"""
Base.occursin(x::Num, y::Num) = occursin(unwrap(x), unwrap(y))
@wrapped function Base.occursin(r::Any, y::Real)
    y = unwrap(y)
    if isequal(r, y)
        return true
    elseif r isa Function
        if r(y)
            return true
        end
    end

    if istree(y)
        return r(operation(y)) ||
                any(y->occursin(r, y), arguments(y))
    else
        return false
    end
end

function filterchildren!(r::Any, y::Union{Num, Symbolic}, acc)
    y = unwrap(y)
    if isequal(r, y)
        push!(acc, y)
        return acc
    elseif r isa Function
        if r(y)
            push!(acc, y)
            return acc
        end
    end

    if istree(y)
        if r(operation(y))
            push!(acc, y)
        end
        foreach(c->filterchildren!(r, c, acc),
                arguments(y))
        return acc
    end
end

filterchildren(r, y) = filterchildren!(r, y, [])

module RewriteHelpers
import Symbolics: is_derivative, filterchildren, unwrap
export replace, occursin, is_derivative,
       filterchildren, unwrap
end
