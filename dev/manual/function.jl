function (ˍ₋out, u)
    #= /home/runner/.julia/packages/SymbolicUtils/ssQsQ/src/code.jl:373 =#
    #= /home/runner/.julia/packages/SymbolicUtils/ssQsQ/src/code.jl:374 =#
    #= /home/runner/.julia/packages/SymbolicUtils/ssQsQ/src/code.jl:375 =#
    begin
        #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:537 =#
        #= /home/runner/.julia/packages/SymbolicUtils/ssQsQ/src/code.jl:422 =# @inbounds begin
                #= /home/runner/.julia/packages/SymbolicUtils/ssQsQ/src/code.jl:418 =#
                ˍ₋out[1] = (+)((getindex)(u, 1), (*)(-1, (getindex)(u, 3)))
                ˍ₋out[2] = (+)((*)(-1, (getindex)(u, 2)), (^)((getindex)(u, 1), 2))
                ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                #= /home/runner/.julia/packages/SymbolicUtils/ssQsQ/src/code.jl:420 =#
                nothing
            end
    end
end