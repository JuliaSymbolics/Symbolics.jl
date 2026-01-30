function (ˍ₋out, u)
    #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:382 =# @inbounds begin
            #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:382 =#
            begin
                #= /home/runner/.julia/packages/SymbolicUtils/kgV0N/src/code.jl:585 =#
                #= /home/runner/.julia/packages/SymbolicUtils/kgV0N/src/code.jl:586 =#
                #= /home/runner/.julia/packages/SymbolicUtils/kgV0N/src/code.jl:587 =#
                #= /home/runner/.julia/packages/SymbolicUtils/kgV0N/src/code.jl:640 =# @inbounds begin
                        #= /home/runner/.julia/packages/SymbolicUtils/kgV0N/src/code.jl:636 =#
                        ˍ₋out[1] = (+)((getindex)(u, 1), (*)(-1, (getindex)(u, 3)))
                        ˍ₋out[2] = (+)((*)(-1, (getindex)(u, 2)), (^)((getindex)(u, 1), 2))
                        ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                        #= /home/runner/.julia/packages/SymbolicUtils/kgV0N/src/code.jl:638 =#
                        ˍ₋out
                    end
            end
        end
end