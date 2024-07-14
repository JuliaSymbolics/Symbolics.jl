const NAMESPACE_SEPARATOR = '₊'

hide_lhs(_) = false

###
### Connection
###
struct Connection
    systems
end
Base.broadcastable(x::Connection) = Ref(x)
Connection() = Connection(nothing)
Base.hash(c::Connection, seed::UInt) = hash(c.systems, (0xc80093537bdc1311 % UInt) ⊻ seed)
hide_lhs(_::Connection) = true

function connect(sys1, sys2, syss...)
    syss = (sys1, sys2, syss...)
    length(unique(nameof, syss)) == length(syss) || error("connect takes distinct systems!")
    Equation(Connection(), Connection(syss)) # the RHS are connected systems
end

function Base.show(io::IO, c::Connection)
    print(io, "connect(")
    if c.systems isa AbstractArray || c.systems isa Tuple
        n = length(c.systems)
        for (i, s) in enumerate(c.systems)
            str = join(split(string(nameof(s)), NAMESPACE_SEPARATOR), '.')
            print(io, str)
            i != n && print(io, ", ")
        end
    end
    print(io, ")")
end

###
### State machine
###
_nameof(s) = nameof(s)
_nameof(s::Union{Int, Symbol}) = s
abstract type StateMachineOperator end
Base.broadcastable(x::StateMachineOperator) = Ref(x)
hide_lhs(_::StateMachineOperator) = true
struct InitialState <: StateMachineOperator
    s
end
Base.show(io::IO, s::InitialState) = print(io, "initial_state(", _nameof(s.s), ")")
initial_state(s) = Equation(InitialState(nothing), InitialState(s))

Base.@kwdef struct Transition{A, B, C} <: StateMachineOperator
    from::A = nothing
    to::B = nothing
    cond::C = nothing
    immediate::Bool = true
    reset::Bool = true
    synchronize::Bool = false
    priority::Int = 1
    function Transition(from, to, cond, immediate, reset, synchronize, priority)
        cond = unwrap(cond)
        new{typeof(from), typeof(to), typeof(cond)}(from, to, cond, immediate,
                                                    reset, synchronize,
                                                    priority)
    end
end
function Base.:(==)(transition1::Transition, transition2::Transition)
    transition1.from == transition2.from &&
    transition1.to == transition2.to &&
    isequal(transition1.cond, transition2.cond) &&
    transition1.immediate == transition2.immediate &&
    transition1.reset == transition2.reset &&
    transition1.synchronize == transition2.synchronize &&
    transition1.priority == transition2.priority
end

"""
    transition(from, to, cond; immediate::Bool = true, reset::Bool = true, synchronize::Bool = false, priority::Int = 1)

Create a transition from state `from` to state `to` that is enabled when transitioncondition `cond` evaluates to `true`.

# Arguments:
- `from`: The source state of the transition.
- `to`: The target state of the transition.
- `cond`: A transition condition that evaluates to a Bool, such as `ticksInState() >= 2`.
- `immediate`: If `true`, the transition will fire at the same tick as it becomes true, if `false`, the actions of the state are evaluated first, and the transition fires during the next tick.
- `reset`: If true, the destination state `to` is reset to its initial condition when the transition fires.
- `synchronize`: If true, the transition will only fire if all sub-state machines in the source state are in their final (terminal) state. A final state is one that has no outgoing transitions.
- `priority`: If a state has more than one outgoing transition, all outgoing transitions must have a unique priority. The transitions are evaluated in priority order, i.e., the transition with priority 1 is evaluated first.
"""
function transition(from, to, cond;
        immediate::Bool = true, reset::Bool = true, synchronize::Bool = false,
        priority::Int = 1)
    Equation(Transition(), Transition(; from, to, cond, immediate, reset,
                                      synchronize, priority))
end
function Base.show(io::IO, s::Transition)
    print(io, _nameof(s.from), " → ", _nameof(s.to), " if (", s.cond, ") [")
    print(io, "immediate: ", Int(s.immediate), ", ")
    print(io, "reset: ", Int(s.reset), ", ")
    print(io, "sync: ", Int(s.synchronize), ", ")
    print(io, "prio: ", s.priority, "]")
end

"""
$(TYPEDEF)

An equality relationship between two expressions.

# Fields
$(FIELDS)
"""
struct Equation
    """The expression on the left-hand side of the equation."""
    lhs
    """The expression on the right-hand side of the equation."""
    rhs
    function Equation(lhs, rhs)
        new(value(lhs), value(rhs))
    end
end
Base.:(==)(a::Equation, b::Equation) = all(isequal.((a.lhs, a.rhs), (b.lhs, b.rhs)))
Base.hash(a::Equation, salt::UInt) = hash(a.lhs, hash(a.rhs, salt))

function Base.show(io::IO, eq::Equation)
    if hide_lhs(eq.lhs)
        show(io, eq.rhs)
    else
        print(io, eq.lhs, " ~ ", eq.rhs)
    end
end

scalarize(eq::Equation) = scalarize(eq.lhs) .~ scalarize(eq.rhs)
SymbolicUtils.simplify(x::Equation; kw...) = simplify(x.lhs; kw...) ~ simplify(x.rhs; kw...)
# ambiguity
for T in [:Pair, :Any]
    @eval function SymbolicUtils.substitute(x::Equation, rules::$T; kw...)
        sub = substituter(rules)
        sub(x.lhs; kw...) ~ sub(x.rhs; kw...)
    end

    @eval function SymbolicUtils.substitute(eqs::Array{Equation}, rules::$T; kw...)
        sub = substituter(rules)
        sub.(lhss(eqs); kw...) .~ sub.(rhss(eqs); kw...)
    end
end

SymbolicUtils.substitute(nums::Array{Num}, rules; kw...) = substituter(rules).(nums; kw...)

lhss(xs) = map(x->x.lhs, xs)
rhss(xs) = map(x->x.rhs, xs)

"""
$(TYPEDSIGNATURES)

Create an [`Equation`](@ref) out of two [`Num`](@ref) instances, or an
`Num` and a `Number`.

# Examples

```jldoctest
julia> using Symbolics

julia> @variables x y;

julia> @variables A[1:3, 1:3] B[1:3, 1:3];

julia> x ~ y
x ~ y

julia> x - y ~ 0
x - y ~ 0

julia> A ~ B
(broadcast(~, A, B))[1:3,1:3]

julia> A .~ 3x
(broadcast(~, A, 3x))[1:3,1:3]
```
"""
function Base.:~(lhs, rhs)
    if (isarraysymbolic(lhs) || isarraysymbolic(rhs)) && ((sl = size(lhs)) != (sr = size(rhs)))
        throw(ArgumentError("Cannot equate an array of different sizes. Got $sl and $sr."))
    else
        Equation(lhs, rhs)
    end
end
for T in [:Num, :Complex, :Number], S in [:Num, :Complex, :Number]
    (T != :Complex && S != :Complex) && continue
    @eval Base.:~(a::$T, b::$S) = let ar = value(real(a)), br = value(real(b)),
                                      ai = value(imag(a)), bi = value(imag(b))
        if ar isa Number && br isa Number && ai isa Number && bi isa Number
            error("Equation $a ~ $b does not contain any symbols")
        elseif ar isa Number && br isa Number
            ai ~ bi
        elseif ai isa Number && bi isa Number
            ar ~ br
        else
            [ar ~ br
            ai ~ bi]
        end
    end
end

canonical_form(eq::Equation) = eq.lhs - eq.rhs ~ 0

get_variables(eq::Equation) = unique(vcat(get_variables(eq.lhs), get_variables(eq.rhs)))

struct ConstrainedEquation
  constraints
  eq
end

function expand_derivatives(eq::Equation, simplify=false)
    return Equation(expand_derivatives(eq.lhs, simplify), expand_derivatives(eq.rhs, simplify))
end
