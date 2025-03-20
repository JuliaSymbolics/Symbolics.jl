function (ˍ₋out, u)
    #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:348 =# @inbounds begin
            #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:348 =#
            begin
                #= /home/runner/.julia/packages/SymbolicUtils/htTbU/src/code.jl:395 =#
                #= /home/runner/.julia/packages/SymbolicUtils/htTbU/src/code.jl:396 =#
                #= /home/runner/.julia/packages/SymbolicUtils/htTbU/src/code.jl:397 =#
                begin
                    #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:564 =#
                    #= /home/runner/.julia/packages/SymbolicUtils/htTbU/src/code.jl:444 =# @inbounds begin
                            #= /home/runner/.julia/packages/SymbolicUtils/htTbU/src/code.jl:440 =#
                            ˍ₋out[1] = (+)((getindex)(u, 1), (*)(-1, (getindex)(u, 3)))
                            ˍ₋out[2] = (+)((*)(-1, (getindex)(u, 2)), (^)((getindex)(u, 1), 2))
                            ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                            #= /home/runner/.julia/packages/SymbolicUtils/htTbU/src/code.jl:442 =#
                            nothing
                        end
                end
            end
        end
end