function (ˍ₋out, u)
    #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:368 =# @inbounds begin
            #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:368 =#
            begin
                #= /home/runner/.julia/packages/SymbolicUtils/a0c2F/src/code.jl:402 =#
                #= /home/runner/.julia/packages/SymbolicUtils/a0c2F/src/code.jl:403 =#
                #= /home/runner/.julia/packages/SymbolicUtils/a0c2F/src/code.jl:404 =#
                #= /home/runner/.julia/packages/SymbolicUtils/a0c2F/src/code.jl:457 =# @inbounds begin
                        #= /home/runner/.julia/packages/SymbolicUtils/a0c2F/src/code.jl:453 =#
                        ˍ₋out[1] = (+)((getindex)(u, 1), (*)(-1, (getindex)(u, 3)))
                        ˍ₋out[2] = (+)((*)(-1, (getindex)(u, 2)), (^)((getindex)(u, 1), 2))
                        ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                        #= /home/runner/.julia/packages/SymbolicUtils/a0c2F/src/code.jl:455 =#
                        ˍ₋out
                    end
            end
        end
end