function (ˍ₋out, u)
    #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:410 =# @inbounds begin
            #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:410 =#
            begin
                #= /home/runner/.julia/packages/SymbolicUtils/Hzfq3/src/code.jl:1219 =# @inbounds begin
                        #= /home/runner/.julia/packages/SymbolicUtils/Hzfq3/src/code.jl:1215 =#
                        ˍ₋out[1] = (+)((getindex)(u, 1), (*)(-1, (getindex)(u, 3)))
                        ˍ₋out[2] = (+)((*)(-1, (getindex)(u, 2)), (^)((getindex)(u, 1), 2))
                        ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                        #= /home/runner/.julia/packages/SymbolicUtils/Hzfq3/src/code.jl:1217 =#
                        ˍ₋out
                    end
            end
        end
end