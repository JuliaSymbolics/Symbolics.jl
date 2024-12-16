function activeState end
function entry end
function ticksInState end
function timeInState end

for (s, T) in [(:timeInState, :Real),
               (:ticksInState, :Integer),
               (:entry, :Bool),
               (:activeState, :Bool)]
    seed = hash(s)
    @eval begin
        $s(x) = wrap(term($s, x))
        SymbolicUtils.promote_symtype(::typeof($s), _...) = $T
        function SymbolicUtils.show_call(io, ::typeof($s), args)
            if isempty(args)
                print(io, $s, "()")
            else
                arg = only(args)
                print(io, $s, "(", arg isa Number ? arg : nameof(arg), ")")
            end
        end
    end
    if s != :activeState
        @eval $s() = wrap(term($s))
    end
end

"""
    timeInState()
    timeInState(state)

Get the time (in seconds) spent in a state in a finite state machine.

When used to query the time spent in the enclosing state, the method without arguments is used, i.e.,
```julia
@mtkmodel FSM begin
    ...
    @equations begin
        var(k+1) ~ timeInState() >= 2 ? 0.0 : var(k)
    end
end
```

If used to query the residence time of another state, the state is passed as an argument.

This operator can be used in both equations and transition conditions.

See also [`ticksInState`](@ref) and [`entry`](@ref)
"""
timeInState

"""
    ticksInState()
    ticksInState(state)

Get the number of ticks spent in a state in a finite state machine.

When used to query the number of ticks spent in the enclosing state, the method without arguments is used, i.e.,
```julia
@mtkmodel FSM begin
    ...
    @equations begin
        var(k+1) ~ ticksInState() >= 2 ? 0.0 : var(k)
    end
end
```

If used to query the number of ticks in another state, the state is passed as an argument.

This operator can be used in both equations and transition conditions.

See also [`timeInState`](@ref) and [`entry`](@ref)
"""
ticksInState

"""
    entry()
    entry(state)

When used in a finite-state machine, this operator returns true at the first tick when the state is active, and false otherwise.

When used to query the entry of the enclosing state, the method without arguments is used, when used to query the entry of another state, the state is passed as an argument.

This can be used to perform a unique action when entering a state.
"""
entry

"""
    activeState(state)

When used in a finite state machine, this operator returns `true` if the queried state is active and false otherwise. 
"""
activeState
