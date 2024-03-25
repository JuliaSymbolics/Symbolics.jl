struct Var
    i::Int
end

"""
Environment: variables which are equivalent as well as
a counter that can be used to create fresh variables
"""
const empty_env = (Base.ImmutableDict{Var, Any}(), 0)
assoc(u, v, s) = Base.ImmutableDict(s, u=>v)

"""
Search: look for Var or return v if not Var
"""
search(s, v::Var) = haskey(s, v) ? search(s, s[v]) : v
search(s, v) = v

"""
Unify:
2 values and a state, return unified state OR
return nothing if unification is not possible
"""
unify(u, v, s) = isequal(u, v) ? s : nothing
unify(u::Var, v::Var, s) = isequal(u, v) ? s : nothing
unify(u::Var, v, s) = assoc(u, v, s)
unify(u, v::Var, s) = assoc(v, u, s)
unify(u::Tuple, v::Tuple, s) = 
    foldl((a, s) -> unify(a..., s),
          init=s, zip(u, v))

"""
Streams:
we use (head, tail) where tail is a tuple.
"""
const empty_stream = ()
struct EmptyList end
struct List
    head
    tail::Union{List, EmptyList}
end
Base.isempty(::EmptyList) = true
Base.isempty(::List) = false
head(l::List) = l.head
tail(l::List) = l.tail

unit(x) = List(x, nothing)

## Goal constuctors
# A goal is a function that takes a state and a counter
# and returns a stream of state counter pairs which unify a and b
function ==ᵤ(a, b)
    function ((st, c))
        st′ = unify(u, v, st)
        !isnothing(s) ? unit((st′, c)) : empty_stream
   end
end

"""
Call the function f with the new variable.
and `f` must return a new goal.
"""
fresh(f) = ((s,c),) -> f(Var(c))((s, c+1))

disj(g1,g2) = sc -> add_streams(g1(sc), g2(sc))
conj(g1,g2) = sc -> bind(g1(sc), g2)

add_streams(s1::Function, s2) = () -> add_streams(s1(), s2)
function add_streams(s1,s2)
    if isempty(s1)
        return s2
    else
        List(s1.head, add_streams(s2[2], s2))
    end
end

mapstream(f, str) = isempty(str) ? () : concat(f(first(str)), mapstream(f, str)

function mapper(stream, goal)
    if isempty(stream)
        return empty_stream
    else
        return add_streams(goal(stream[1]), bind(stream[2], goal))
    end
end
bind(stream::Function, goal)  = () -> bind(stream(), goal)
