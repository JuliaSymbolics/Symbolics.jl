function (ˍ₋out, u)
    #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:368 =# @inbounds begin
            #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:368 =#
            begin
                #= /home/runner/.julia/packages/SymbolicUtils/KmZ71/src/code.jl:409 =#
                #= /home/runner/.julia/packages/SymbolicUtils/KmZ71/src/code.jl:410 =#
                #= /home/runner/.julia/packages/SymbolicUtils/KmZ71/src/code.jl:411 =#
                #= /home/runner/.julia/packages/SymbolicUtils/KmZ71/src/code.jl:464 =# @inbounds begin
                        #= /home/runner/.julia/packages/SymbolicUtils/KmZ71/src/code.jl:460 =#
                        ˍ₋out[1] = (+)((getindex)(u, 1), (*)(-1, (getindex)(u, 3)))
                        ˍ₋out[2] = (+)((*)(-1, (getindex)(u, 2)), (^)((getindex)(u, 1), 2))
                        ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                        #= /home/runner/.julia/packages/SymbolicUtils/KmZ71/src/code.jl:462 =#
                        ˍ₋out
                    end
            end
        end
end