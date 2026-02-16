function (ˍ₋out, u)
    #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:383 =# @inbounds begin
            #= /home/runner/work/Symbolics.jl/Symbolics.jl/src/build_function.jl:383 =#
            begin
                #= /home/runner/.julia/packages/SymbolicUtils/DAriZ/src/code.jl:711 =#
                #= /home/runner/.julia/packages/SymbolicUtils/DAriZ/src/code.jl:712 =#
                #= /home/runner/.julia/packages/SymbolicUtils/DAriZ/src/code.jl:713 =#
                #= /home/runner/.julia/packages/SymbolicUtils/DAriZ/src/code.jl:813 =# @inbounds begin
                        #= /home/runner/.julia/packages/SymbolicUtils/DAriZ/src/code.jl:809 =#
                        ˍ₋out[1] = (+)((getindex)(u, 1), (*)(-1, (getindex)(u, 3)))
                        ˍ₋out[2] = (+)((*)(-1, (getindex)(u, 2)), (^)((getindex)(u, 1), 2))
                        ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                        #= /home/runner/.julia/packages/SymbolicUtils/DAriZ/src/code.jl:811 =#
                        ˍ₋out
                    end
            end
        end
end