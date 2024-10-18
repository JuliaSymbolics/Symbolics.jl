function (ˍ₋out, u)
    #= /home/runner/.julia/packages/SymbolicUtils/jf8aQ/src/code.jl:385 =#
    #= /home/runner/.julia/packages/SymbolicUtils/jf8aQ/src/code.jl:386 =#
    #= /home/runner/.julia/packages/SymbolicUtils/jf8aQ/src/code.jl:387 =#
    begin
        #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:546 =#
        #= /home/runner/.julia/packages/SymbolicUtils/jf8aQ/src/code.jl:434 =# @inbounds begin
                #= /home/runner/.julia/packages/SymbolicUtils/jf8aQ/src/code.jl:430 =#
                ˍ₋out[1] = (+)((getindex)(u, 1), (*)(-1, (getindex)(u, 3)))
                ˍ₋out[2] = (+)((*)(-1, (getindex)(u, 2)), (^)((getindex)(u, 1), 2))
                ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                #= /home/runner/.julia/packages/SymbolicUtils/jf8aQ/src/code.jl:432 =#
                nothing
            end
    end
end