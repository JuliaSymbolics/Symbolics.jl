function (ˍ₋out, u)
    #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:383 =# @inbounds begin
            #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:383 =#
            begin
                #= /home/runner/.julia/packages/SymbolicUtils/u8OUn/src/code.jl:1148 =# @inbounds begin
                        #= /home/runner/.julia/packages/SymbolicUtils/u8OUn/src/code.jl:1144 =#
                        ˍ₋out[1] = (+)((getindex)(u, 1), (*)(-1, (getindex)(u, 3)))
                        ˍ₋out[2] = (+)((*)(-1, (getindex)(u, 2)), (^)((getindex)(u, 1), 2))
                        ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                        #= /home/runner/.julia/packages/SymbolicUtils/u8OUn/src/code.jl:1146 =#
                        ˍ₋out
                    end
            end
        end
end