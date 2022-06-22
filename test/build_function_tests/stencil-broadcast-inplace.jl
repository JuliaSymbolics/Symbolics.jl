:(function (x,)
      let ˍ₋out = zeros(Float64, map(length, (1:6, 1:6)))
          begin
              ˍ₋out_2_input_1 = (broadcast)(+, x, (adjoint)(x))
              ˍ₋out_1 = (view)(ˍ₋out, 1:6, 1:6)
              ˍ₋out_1 .= 0
              ˍ₋out_2 = (view)(ˍ₋out, 2:5, 2:5)
              for (j, j′) = zip(1:4, reset_to_one(1:4))
                  for (i, i′) = zip(1:4, reset_to_one(1:4))
                      ˍ₋out_2[i′, j′] = (+)(ˍ₋out_2[i′, j′], (+)(1, (getindex)(ˍ₋out_2_input_1, i, j)))
                  end
              end
          end
          ˍ₋out
      end
  end)