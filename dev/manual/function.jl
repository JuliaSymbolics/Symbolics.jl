function (ˍ₋out, u)
    #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:347 =# @inbounds begin
            #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:347 =#
            begin
                #= /home/runner/.julia/packages/SymbolicUtils/2cE5G/src/code.jl:388 =#
                #= /home/runner/.julia/packages/SymbolicUtils/2cE5G/src/code.jl:389 =#
                #= /home/runner/.julia/packages/SymbolicUtils/2cE5G/src/code.jl:390 =#
                begin
                    #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:563 =#
                    #= /home/runner/.julia/packages/SymbolicUtils/2cE5G/src/code.jl:437 =# @inbounds begin
                            #= /home/runner/.julia/packages/SymbolicUtils/2cE5G/src/code.jl:433 =#
                            ˍ₋out[1] = (+)((getindex)(u, 1), (*)(-1, (getindex)(u, 3)))
                            ˍ₋out[2] = (+)((*)(-1, (getindex)(u, 2)), (^)((getindex)(u, 1), 2))
                            ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                            #= /home/runner/.julia/packages/SymbolicUtils/2cE5G/src/code.jl:435 =#
                            nothing
                        end
                end
            end
        end
end