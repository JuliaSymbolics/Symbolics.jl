function (ˍ₋out, u)
    #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:383 =# @inbounds begin
            #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:383 =#
            begin
                #= /home/runner/.julia/packages/SymbolicUtils/vtVwT/src/code.jl:584 =#
                #= /home/runner/.julia/packages/SymbolicUtils/vtVwT/src/code.jl:585 =#
                #= /home/runner/.julia/packages/SymbolicUtils/vtVwT/src/code.jl:586 =#
                #= /home/runner/.julia/packages/SymbolicUtils/vtVwT/src/code.jl:639 =# @inbounds begin
                        #= /home/runner/.julia/packages/SymbolicUtils/vtVwT/src/code.jl:635 =#
                        ˍ₋out[1] = (+)((getindex)(u, 1), (*)(-1, (getindex)(u, 3)))
                        ˍ₋out[2] = (+)((*)(-1, (getindex)(u, 2)), (^)((getindex)(u, 1), 2))
                        ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                        #= /home/runner/.julia/packages/SymbolicUtils/vtVwT/src/code.jl:637 =#
                        ˍ₋out
                    end
            end
        end
end