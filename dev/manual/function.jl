function (ˍ₋out, u)
    #= /home/runner/.julia/packages/SymbolicUtils/1JRDc/src/code.jl:350 =#
    #= /home/runner/.julia/packages/SymbolicUtils/1JRDc/src/code.jl:351 =#
    #= /home/runner/.julia/packages/SymbolicUtils/1JRDc/src/code.jl:352 =#
    begin
        #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:520 =#
        #= /home/runner/.julia/packages/SymbolicUtils/1JRDc/src/code.jl:399 =# @inbounds begin
                #= /home/runner/.julia/packages/SymbolicUtils/1JRDc/src/code.jl:395 =#
                ˍ₋out[1] = (+)((*)(-1, (getindex)(u, 3)), (getindex)(u, 1))
                ˍ₋out[2] = (+)((^)((getindex)(u, 1), 2), (*)(-1, (getindex)(u, 2)))
                ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                #= /home/runner/.julia/packages/SymbolicUtils/1JRDc/src/code.jl:397 =#
                nothing
            end
    end
end