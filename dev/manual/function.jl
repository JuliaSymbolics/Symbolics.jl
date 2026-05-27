function (ˍ₋out, u)
    #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:383 =# @inbounds begin
            #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:383 =#
            begin
                #= /home/runner/.julia/packages/SymbolicUtils/bdo1i/src/code.jl:1118 =# @inbounds begin
                        #= /home/runner/.julia/packages/SymbolicUtils/bdo1i/src/code.jl:1114 =#
                        ˍ₋out[1] = (+)((getindex)(u, 1), (*)(-1, (getindex)(u, 3)))
                        ˍ₋out[2] = (+)((*)(-1, (getindex)(u, 2)), (^)((getindex)(u, 1), 2))
                        ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                        #= /home/runner/.julia/packages/SymbolicUtils/bdo1i/src/code.jl:1116 =#
                        ˍ₋out
                    end
            end
        end
end