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
