function (ˍ₋out, u)
    #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:348 =# @inbounds begin
            #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:348 =#
            begin
                #= /home/runner/.julia/packages/SymbolicUtils/6fncq/src/code.jl:389 =#
                #= /home/runner/.julia/packages/SymbolicUtils/6fncq/src/code.jl:390 =#
                #= /home/runner/.julia/packages/SymbolicUtils/6fncq/src/code.jl:391 =#
                begin
                    #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:564 =#
                    #= /home/runner/.julia/packages/SymbolicUtils/6fncq/src/code.jl:438 =# @inbounds begin
                            #= /home/runner/.julia/packages/SymbolicUtils/6fncq/src/code.jl:434 =#
                            ˍ₋out[1] = (+)((getindex)(u, 1), (*)(-1, (getindex)(u, 3)))
                            ˍ₋out[2] = (+)((*)(-1, (getindex)(u, 2)), (^)((getindex)(u, 1), 2))
                            ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                            #= /home/runner/.julia/packages/SymbolicUtils/6fncq/src/code.jl:436 =#
                            nothing
                        end
                end
            end
        end
end