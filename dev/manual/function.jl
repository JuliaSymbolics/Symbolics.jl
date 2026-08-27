function (ˍ₋out, u)
    #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:410 =# @inbounds begin
            #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:410 =#
            begin
                #= /home/runner/.julia/packages/SymbolicUtils/xkUyd/src/code.jl:1264 =# @inbounds begin
                        #= /home/runner/.julia/packages/SymbolicUtils/xkUyd/src/code.jl:1260 =#
                        ˍ₋out[1] = (+)((getindex)(u, 1), (*)(-1, (getindex)(u, 3)))
                        ˍ₋out[2] = (+)((*)(-1, (getindex)(u, 2)), (^)((getindex)(u, 1), 2))
                        ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                        #= /home/runner/.julia/packages/SymbolicUtils/xkUyd/src/code.jl:1262 =#
                        ˍ₋out
                    end
            end
        end
end