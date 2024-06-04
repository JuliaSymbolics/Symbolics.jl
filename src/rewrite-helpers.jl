"""
replacenode(expr::Symbolic, rules...)
Walk the expression and replacenode subexpressions according to `rules`. `rules`
could be rules constructed with `@rule`, a function, or a pair where the
left hand side is matched with equality (using `isequal`) and is replacenoded by the right hand side.

Rules will be applied left-to-right simultaneously,
so only one pattern will be applied to any subexpression,
and the patterns will only be applied to the input text,
not the replacenodements.

Set `fixpoint = true` to repeatedly apply rules until no
change to the expression remains to be made.
"""
function replacenode(expr::Num, r::Pair, rules::Pair...; fixpoint = false)
    _replacenode(unwrap(expr), r, rules...)
end
# Fix ambiguity
replacenode(expr::Num, rules...; fixpoint = false) = _replacenode(unwrap(expr), rules...; fixpoint)
replacenode(expr::Symbolic, rules...; fixpoint = false) = _replacenode(unwrap(expr), rules...; fixpoint)
replacenode(expr::Symbolic, r::Pair, rules::Pair...; fixpoint = false) = _replacenode(expr, r, rules...; fixpoint)
replacenode(expr::Number, rules...; fixpoint = false) = expr
replacenode(expr::Number, r::Pair, rules::Pair...; fixpoint = false) = expr

function _replacenode(expr::Symbolic, rules...; fixpoint = false)
    rs = map(r -> r isa Pair ? (x -> isequal(x, r[1]) ? r[2] : nothing) : r, rules)
    R = Prewalk(Chain(rs))
    if fixpoint
        Fixpoint(R)(expr)
    else
        R(expr)
    end
end

"""
    hasnode(c, x)
Returns true if any part of `x` fufills the condition given in c. c can be a function or an expression.
If it is a function, returns true if x is true for any part of x. If c is an expression, returns
true if x contains c.

Examples:
```julia
@syms x y
hasnode(x, log(x) + x + 1) # returns `true`.
hasnode(x, log(y) + y + 1) # returns `false`.
```

```julia
@variables t X(t)
D = Differential(t)
hasnode(Symbolics.is_derivative, X + D(X) + D(X^2)) # returns `true`.
```
"""
function hasnode(r::Function, y::Union{Num, Symbolic})
    _hasnode(r, y)
end
hasnode(r::Num, y::Union{Num, Symbolic}) = occursin(unwrap(r), unwrap(y))
hasnode(r::Symbolic, y::Union{Num, Symbolic}) = occursin(unwrap(r), unwrap(y))
hasnode(r::Union{Num, Symbolic, Function}, y::Number) = false

function _hasnode(r, y)
    y = unwrap(y)
    if r isa Function
        if r(y)
            return true
        end
    end

    if iscall(y)
        return r(operation(y)) ||
                any(y->_hasnode(r, y), arguments(y))
    else
        return false
    end
end

"""
filterchildren(c, x)
Returns all parts of `x` that fufills the condition given in c. c can be a function or an expression.
If it is a function, returns everything for which the function is `true`. If c is an expression, returns
all expressions that matches it.

Examples:
```julia
@syms x
Symbolics.filterchildren(x, log(x) + x + 1)
```
returns `[x, x]`

```julia
@variables t X(t)
D = Differential(t)
Symbolics.filterchildren(Symbolics.is_derivative, X + D(X) + D(X^2))
```
returns `[Differential(t)(X(t)^2), Differential(t)(X(t))]`
"""
filterchildren(r, y) = filterchildren!(r, y, [])

function filterchildren!(r::Any, y, acc)
    y = unwrap(y)
    r = unwrap(r)
    if isequal(r, y)
        push!(acc, y)
        return acc
    elseif r isa Function
        if r(y)
            push!(acc, y)
            return acc
        end
    end

    if iscall(y)
        if isequal(r, operation(y))
            push!(acc, operation(y))
        elseif r isa Function && r(operation(y))
            push!(acc, operation(y))
        end
        foreach(c->filterchildren!(r, c, acc),
                arguments(y))
        return acc
    end
end

module RewriteHelpers
import Symbolics: replacenode, hasnode, filterchildren, unwrap
export replacenode, hasnode, filterchildren, unwrap

end
