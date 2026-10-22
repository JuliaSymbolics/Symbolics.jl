function (ˍ₋out, u)
    #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:456 =# @inbounds begin
            #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:456 =#
            begin
                #= /home/runner/.julia/packages/SymbolicUtils/HQQQT/src/code.jl:1298 =# @inbounds begin
                        #= /home/runner/.julia/packages/SymbolicUtils/HQQQT/src/code.jl:1294 =#
                        ˍ₋out[1] = (+)((getindex)(u, 1), (*)(-1, (getindex)(u, 3)))
                        ˍ₋out[2] = (+)((*)(-1, (getindex)(u, 2)), (Base.literal_pow)(^, (getindex)(u, 1), (Val)(2)))
                        ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                        #= /home/runner/.julia/packages/SymbolicUtils/HQQQT/src/code.jl:1296 =#
                        ˍ₋out
                    end
            end
        end
end