function (ˍ₋out, u)
    #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:465 =# @inbounds begin
            #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:465 =#
            begin
                #= /home/runner/.julia/packages/SymbolicUtils/S0V4b/src/code.jl:1312 =# @inbounds begin
                        #= /home/runner/.julia/packages/SymbolicUtils/S0V4b/src/code.jl:1308 =#
                        ˍ₋out[1] = (+)((getindex)(u, 1), (*)(-1, (getindex)(u, 3)))
                        ˍ₋out[2] = (+)((*)(-1, (getindex)(u, 2)), (Base.literal_pow)(^, (getindex)(u, 1), (Val)(2)))
                        ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                        #= /home/runner/.julia/packages/SymbolicUtils/S0V4b/src/code.jl:1310 =#
                        ˍ₋out
                    end
            end
        end
end