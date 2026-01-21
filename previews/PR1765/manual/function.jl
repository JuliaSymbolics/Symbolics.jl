function (ˍ₋out, u)
    #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:382 =# @inbounds begin
            #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:382 =#
            begin
                #= /home/runner/.julia/packages/SymbolicUtils/WcZJi/src/code.jl:540 =#
                #= /home/runner/.julia/packages/SymbolicUtils/WcZJi/src/code.jl:541 =#
                #= /home/runner/.julia/packages/SymbolicUtils/WcZJi/src/code.jl:542 =#
                #= /home/runner/.julia/packages/SymbolicUtils/WcZJi/src/code.jl:595 =# @inbounds begin
                        #= /home/runner/.julia/packages/SymbolicUtils/WcZJi/src/code.jl:591 =#
                        ˍ₋out[1] = (+)((getindex)(u, 1), (*)(-1, (getindex)(u, 3)))
                        ˍ₋out[2] = (+)((*)(-1, (getindex)(u, 2)), (^)((getindex)(u, 1), 2))
                        ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                        #= /home/runner/.julia/packages/SymbolicUtils/WcZJi/src/code.jl:593 =#
                        ˍ₋out
                    end
            end
        end
end