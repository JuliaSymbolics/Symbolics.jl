function (ˍ₋out, u)
    #= /home/runner/.julia/packages/SymbolicUtils/qulQp/src/code.jl:349 =#
    #= /home/runner/.julia/packages/SymbolicUtils/qulQp/src/code.jl:350 =#
    #= /home/runner/.julia/packages/SymbolicUtils/qulQp/src/code.jl:351 =#
    begin
        #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:519 =#
        #= /home/runner/.julia/packages/SymbolicUtils/qulQp/src/code.jl:398 =# @inbounds begin
                #= /home/runner/.julia/packages/SymbolicUtils/qulQp/src/code.jl:394 =#
                ˍ₋out[1] = (+)((*)(-1, (getindex)(u, 3)), (getindex)(u, 1))
                ˍ₋out[2] = (+)((^)((getindex)(u, 1), 2), (*)(-1, (getindex)(u, 2)))
                ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                #= /home/runner/.julia/packages/SymbolicUtils/qulQp/src/code.jl:396 =#
                nothing
            end
    end
end